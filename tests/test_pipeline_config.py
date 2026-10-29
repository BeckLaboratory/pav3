"""Tests for locating and reading the PAV configuration file."""
from __future__ import annotations

import pytest

from pav3.params import PavParams
from pav3.pipeline import get_pav_asm_table_path, get_pav_config_path, read_pav_config

CONFIG_JSON = """{
  "reference": "/path/to/hg38.no_alt.fa.gz",
  "min_mapq": 10,
  "vcf_cpx_symbolic": true,
  "merge_params": {"insdel": {"ro_min": 0.5}}
}
"""

CONFIG_YAML = """reference: /path/to/hg38.no_alt.fa.gz
min_mapq: 10
vcf_cpx_symbolic: true
merge_params:
  insdel:
    ro_min: 0.5
"""


def test_config_path_yaml(tmp_path):
    """A YAML config is found."""
    (tmp_path / 'pav.yaml').write_text(CONFIG_YAML)

    assert get_pav_config_path(root_path=tmp_path) == tmp_path / 'pav.yaml'


def test_config_path_yml(tmp_path):
    """A YAML config with the ".yml" extension is found."""
    (tmp_path / 'pav.yml').write_text(CONFIG_YAML)

    assert get_pav_config_path(root_path=tmp_path) == tmp_path / 'pav.yml'


def test_config_path_json_fallback(tmp_path):
    """Without a YAML config, the JSON config is used."""
    (tmp_path / 'pav.json').write_text(CONFIG_JSON)

    assert get_pav_config_path(root_path=tmp_path) == tmp_path / 'pav.json'


def test_config_path_yaml_and_json(tmp_path, capsys):
    """With both YAML and JSON configs, the YAML config is used and a warning names the ignored JSON config."""
    (tmp_path / 'pav.yaml').write_text(CONFIG_YAML)
    (tmp_path / 'pav.json').write_text(CONFIG_JSON)

    assert get_pav_config_path(root_path=tmp_path) == tmp_path / 'pav.yaml'

    out = capsys.readouterr().out
    assert 'WARNING' in out
    assert str(tmp_path / 'pav.json') in out


def test_config_path_missing(tmp_path):
    """A missing config raises FileNotFoundError naming the expected config."""
    with pytest.raises(FileNotFoundError, match='pav.yaml'):
        get_pav_config_path(root_path=tmp_path)


def test_read_config_yaml_matches_json(tmp_path):
    """YAML and JSON configs with the same values produce the same config and parameters."""
    yaml_path = tmp_path / 'pav.yaml'
    json_path = tmp_path / 'pav.json'

    yaml_path.write_text(CONFIG_YAML)
    json_path.write_text(CONFIG_JSON)

    config_yaml = read_pav_config(yaml_path)
    config_json = read_pav_config(json_path)

    assert config_yaml == config_json

    params = PavParams('asm', config_yaml)
    assert params.min_mapq == 10
    assert params.vcf_cpx_symbolic is True
    assert params.merge_params == {'insdel': {'ro_min': 0.5}}


def test_read_config_json_txt(tmp_path):
    """A ".json.txt" config is parsed as JSON."""
    path = tmp_path / 'pav.json.txt'
    path.write_text(CONFIG_JSON)

    assert read_pav_config(path)['min_mapq'] == 10


@pytest.mark.parametrize(
    'filename, content',
    [
        ('pav.yaml', 'reference: [unclosed\n'),
        ('pav.json', '{"reference": }'),
        ('pav.yaml', '- reference\n'),
        ('pav.yaml', ''),
    ]
)
def test_read_config_invalid(tmp_path, filename, content):
    """Malformed configs and configs without a top-level mapping raise ValueError."""
    path = tmp_path / filename
    path.write_text(content)

    with pytest.raises(ValueError, match=filename):
        read_pav_config(path)


@pytest.mark.parametrize(
    'present, expected',
    [
        (['assemblies.csv', 'assemblies.xlsx', 'assemblies.tsv'], 'assemblies.csv'),
        (['assemblies.xlsx', 'assemblies.tsv'], 'assemblies.xlsx'),
        (['assemblies.tsv'], 'assemblies.tsv'),
    ]
)
def test_asm_table_path_order(tmp_path, present, expected):
    """Default assembly tables are searched in order: CSV, Excel, then TSV."""
    for filename in present:
        (tmp_path / filename).write_text('')

    assert get_pav_asm_table_path(root_path=tmp_path) == tmp_path / expected


def test_asm_table_path_config(tmp_path):
    """An assembly table set in the configuration is used over the defaults."""
    (tmp_path / 'assemblies.csv').write_text('')
    (tmp_path / 'samples.parquet').write_text('')

    assert get_pav_asm_table_path({'assembly_table': 'samples.parquet'}, tmp_path) == tmp_path / 'samples.parquet'


def test_asm_table_path_missing(tmp_path):
    """Without an assembly table, the error names each default path."""
    with pytest.raises(FileNotFoundError, match='assemblies.csv.*assemblies.xlsx.*assemblies.tsv'):
        get_pav_asm_table_path(root_path=tmp_path)
