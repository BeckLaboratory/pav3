"""Tests for rule performance logging (config parameter "perf_log")."""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import polars as pl
import pytest

from pav3 import const
from pav3.params import PavParams
from pav3.workflow import perf

META_COLS = ['rule', 'asm_name', 'hap', 'wildcards']


@pytest.mark.parametrize(
    'val, expected',
    [
        (None, 0),
        (True, const.DEFAULT_PERF_LOG_INTERVAL),
        (False, 0),
        (0, 0),
        (1, 1),
        (45, 45),
        ('true', const.DEFAULT_PERF_LOG_INTERVAL),
        ('False', 0),
        ('1', 1),  # Integer before boolean: 1 second, not True
        ('0', 0),
        (' 60 ', 60),
    ]
)
def test_perf_log_interval(val, expected):
    """Booleans select the default interval or disable logging, integers set the interval."""
    config = {} if val is None else {'perf_log': val}

    assert PavParams('asm', config).perf_log == expected


@pytest.mark.parametrize('val', [-1, '-5', 1.5, '1.5', 'maybe', [30]])
def test_perf_log_interval_illegal(val):
    """Negative, fractional, and unrecognized values are errors."""
    with pytest.raises(ValueError, match='perf_log'):
        _ = PavParams('asm', {'perf_log': val}).perf_log


def test_perf_log_interval_asm_override():
    """The assembly table config column overrides the config file for that assembly."""
    asm_table = pl.DataFrame({'name': ['asm1', 'asm2'], 'config': ['perf_log=10', None]})

    assert PavParams('asm1', {'perf_log': True}, asm_table).perf_log == 10
    assert PavParams('asm2', {'perf_log': True}, asm_table).perf_log == const.DEFAULT_PERF_LOG_INTERVAL


@pytest.mark.parametrize(
    'rule, wildcards, expected',
    [
        (
            'call_intra_snv_insdel', {'asm_name': 'HG00171', 'hap': 'h1'},
            'log/perf/HG00171/call_intra_snv_insdel_h1.parquet'
        ),
        ('call_tables', {'asm_name': 'HG00171', 'vartype': 'insdel'}, 'log/perf/HG00171/call_tables_insdel.parquet'),
        ('call_vcf', {'asm_name': 'HG00171'}, 'log/perf/HG00171/call_vcf.parquet'),
        ('data_ref_fofn', {}, 'log/perf/data_ref_fofn.parquet'),
    ]
)
def test_perf_log_path(rule, wildcards, expected):
    """Wildcards other than asm_name are appended to the rule name, and asm_name is the directory."""
    assert perf.perf_log_path(rule, wildcards) == Path(expected)


def test_perf_log_disabled(tmp_path, monkeypatch):
    """Disabled logging writes nothing."""
    monkeypatch.chdir(tmp_path)

    with perf.perf_log('call_vcf', {'asm_name': 'asm'}, {'perf_log': False}):
        pass

    assert not (tmp_path / 'log').exists()


def test_perf_log_enabled(tmp_path, monkeypatch):
    """A rule's table identifies the rule, measures from the start of the rule, and ends with a final sample."""
    monkeypatch.chdir(tmp_path)

    with perf.perf_log('call_inter', {'asm_name': 'asm', 'hap': 'h1'}, {'perf_log': 1}) as perf_obj:
        assert isinstance(perf_obj, perf.PerfLog)
        time.sleep(0.2)

    df = pl.read_parquet(tmp_path / 'log/perf/asm/call_inter_h1.parquet')

    assert df.columns == META_COLS + list(perf.build_schema(perf.parse_columns(perf.PERF_COLUMNS)))
    assert df.select(META_COLS).unique().rows() == [('call_inter', 'asm', 'h1', 'hap=h1')]

    # First sample on entry and a final sample on exit
    assert df.height == 2
    assert df['elapsed_s'][0] < 0.1
    assert df['elapsed_s'][-1] >= 0.2

    # Cumulative counters start at the rule, not at the start of the process
    assert df['cpu_s'][0] == 0.0
    assert df['cores_n'][0] == 0.0


def test_perf_log_samples_tree(tmp_path):
    """Samples follow the schedule, include child processes, and keep the CPU of children that exited."""
    path = tmp_path / 'perf.parquet'
    cmd = [sys.executable, '-c', 'import time\nt = time.time()\nwhile time.time() - t < 0.6: pass']

    with perf.PerfLog(path, 0.1):
        subprocess.run(cmd, check=True)

    df = pl.read_parquet(path)

    assert df.height >= 5
    assert df['proc_n'].max() >= 2
    assert df['proc_n'][-1] == 1
    assert df['cpu_s'][-1] >= 0.4
    assert df['cpu_s'].diff().drop_nulls().min() >= 0.0


def test_perf_log_written_while_running(tmp_path):
    """The table is rewritten after each sample, so a job killed before exiting keeps its samples."""
    path = tmp_path / 'perf.parquet'

    with perf.PerfLog(path, 0.1):
        assert pl.read_parquet(path).height == 1

        time.sleep(0.35)

        assert pl.read_parquet(path).height >= 3

    assert not list(tmp_path.glob('.*.tmp'))


def test_perf_log_exception(tmp_path):
    """The table is written when the block fails, and the block's exception is raised."""
    path = tmp_path / 'perf.parquet'

    with pytest.raises(RuntimeError, match='rule failed'):
        with perf.PerfLog(path, 1, meta={'rule': 'test'}):
            raise RuntimeError('rule failed')

    df = pl.read_parquet(path)

    assert df.height == 2
    assert df['rule'].to_list() == ['test', 'test']


def test_perf_log_unwritable(tmp_path):
    """A table that cannot be written fails on entry before the block runs."""
    (tmp_path / 'file').touch()
    ran = False

    with pytest.raises(OSError):
        with perf.PerfLog(tmp_path / 'file' / 'perf.parquet', 1):
            ran = True

    assert not ran


def test_perf_log_args():
    """A sampling interval must be positive, and a meta column cannot replace a sampled column."""
    with pytest.raises(ValueError, match='positive'):
        perf.PerfLog('perf.parquet', 0)

    with pytest.raises(ValueError, match='collides'):
        perf.PerfLog('perf.parquet', 1, meta={'rss_gb': 'x'})
