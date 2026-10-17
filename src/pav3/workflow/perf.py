"""Sample the resource usage of a workflow rule while it runs.

Adapted from the resource sampler shared by the variant_discovery pipelines (``<tool>_perf.py``, itself adapted from
becktools' ``run-with-perf.py``). Those pipelines launch a command under a sampler process. A PAV rule is Python code
running in its own Snakemake job process, so here the sampler is a background thread that follows that process and
every descendant it still owns (e.g. minimap2 launched through ``shell()``). The column definitions, the process-tree
accounting (:class:`TreeMonitor`), and the sampling schedule are carried over unchanged, so the measurement columns
match the variant_discovery tables. The command runner, signal relay, command line, and output formats other than
parquet are not carried over.

Each rule from input assemblies to the VCF runs inside :func:`perf_log`. When config parameter ``perf_log`` enables
it, the rule writes a table to ``log/perf/{asm_name}/{rule}_{wildcards}.parquet`` with constant columns ``rule``,
``asm_name``, ``hap``, and ``wildcards`` ahead of the measurements, so the tables identify themselves once they are
concatenated (``pl.read_parquet('log/perf/**/*.parquet')``).

Three behaviors differ from the variant_discovery sampler:

* The table is rewritten (atomically) after every sample instead of once when the command exits. A job the scheduler
  kills (SIGKILL for exceeding its memory) takes the sampler thread with it, and the job that ran out of memory is
  the one whose numbers are wanted most. A killed job keeps every sample up to the last interval.
* A final sample is taken when the rule ends. A sampler process cannot sample a command that has exited, but a thread
  can sample its own rule, so the last row reports the rule's full run time and CPU seconds.
* Cumulative counters (``cpu_s``, ``read_*``, ``write_*``) start at 0 when the rule starts. A launched command starts
  near 0 on its own, but the Snakemake job process has already spent CPU and I/O starting Snakemake.

Following the variant_discovery convention, a column name carries its unit as a suffix: ``_b`` bytes, ``_mb``
mebibytes, ``_gb`` gibibytes, ``_s`` seconds, ``_n`` a count or a fraction of a resource where ``1.0`` is 100% of one
unit of it (so ``cores_n`` of ``2.5`` means two and a half cores' worth of CPU, not 2.5%), and ``_bps`` bytes per
second. PAV writes every column (:const:`PERF_COLUMNS`) with byte values in GiB.

Memory is summed over the process tree. ``pss_*`` (proportional set size, Linux only) divides each shared page among
the processes mapping it, so summing it over the tree neither double-counts memory shared *within* the tree nor
charges the tree for the whole of a page it shares with an unrelated process; it is the most faithful measure of what
the rule actually costs the machine. ``rss_*`` double-counts pages shared between processes in the tree, and
``vss_*`` counts address space that may never be resident. On platforms without PSS the ``pss_*`` and ``swap_*``
columns are null.

CPU is reported as ``cores_n``: the CPU seconds consumed by the tree between two samples divided by the wall-clock
seconds between them. A process that exits does not take its CPU time out of the total with it. Two independent lower
bounds are kept and the larger is used: the time of descendants each live process has already reaped
(``children_user`` / ``children_system``, which also catches children that came and went entirely between two
samples, but which only Linux reports), and the time of every process a sample has actually seen alive. The counter
is held monotonic, so a process sitting unreaped between samples can never produce a negative rate.

The tree is rooted at the Snakemake job process. Snakemake runs each job with a ``run`` block in its own process, so
the tree is the rule and nothing else, unless Snakemake is run with ``--force-use-threads``, which runs jobs as
threads of the main process and puts every concurrent job in the same tree.

Limitations inherent to sampling a process tree: a descendant that is reparented away (daemonized) leaves the tree and
stops being counted, and a spike shorter than the sampling interval can be missed entirely. Only the CPU accounting is
immune to the latter.
"""

__all__ = [
    'PERF_COLUMNS',
    'PERF_LOG_DIR',
    'PerfLog',
    'TreeMonitor',
    'parse_columns',
    'perf_log',
    'perf_log_path',
]

import contextlib
import logging
import os
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from types import TracebackType
from typing import Any, Optional

import polars as pl
import psutil

from ..params import PavParams

logger = logging.getLogger(__name__)

PERF_LOG_DIR: Path = Path('log/perf')
"""Directory performance tables are written to, relative to the PAV run directory."""

PERF_COLUMNS: str = 'all_gb'
"""Columns written by PAV rules (every column, byte values in GiB), the same as the variant_discovery pipelines."""

#: Fields making up the default columns, in output order. These are always present: ``columns`` adds to them rather
#: than replacing them. Byte fields among them are rendered in :const:`DEFAULT_UNIT` unless a ``default_*`` or
#: ``all_*`` keyword says otherwise.
DEFAULT_FIELDS: tuple[str, ...] = ("time", "elapsed_s", "pss", "rss", "vss", "cores_n", "proc_n")

#: Every field in canonical output order, used by the ``all`` keyword. The default fields lead,
#: so adding ``all`` appends columns rather than reordering the ones already there.
ALL_FIELDS: tuple[str, ...] = DEFAULT_FIELDS + (
    "cpu_s",
    "uss",
    "swap",
    "threads_n",
    "fd_n",
    "mem_n",
    "read",
    "write",
    "read_bps",
    "write_bps",
)

#: Byte unit used when nothing more specific selects one.
DEFAULT_UNIT: str = "b"

#: Keywords naming a group of columns rather than a single one. Each accepts a byte-unit
#: suffix that sets the unit for every byte column in its group.
GROUP_KEYWORDS: tuple[str, ...] = ("all", "default")

#: Sample fields holding a byte count. Each is requested with a byte-unit suffix from
#: :const:`BYTE_UNITS` (e.g. ``rss_b``, ``rss_mb``, ``rss_gb``) and has no bare form.
BYTE_FIELDS: dict[str, str] = {
    "pss": "Proportional set size: shared pages divided among the processes mapping them "
           "(Linux only; null elsewhere). The most accurate measure of the tree's memory.",
    "rss": "Resident set size: physical memory in use, double-counting pages shared between "
           "processes in the tree.",
    "vss": "Virtual set size: total address space mapped, resident or not.",
    "uss": "Unique set size: memory private to the tree, i.e. what would be freed if it exited.",
    "swap": "Swapped-out memory (Linux only; null elsewhere).",
    "read": "Cumulative bytes read by the tree since the rule started (Linux only; null elsewhere).",
    "write": "Cumulative bytes written by the tree since the rule started (Linux only; null "
             "elsewhere).",
}

#: Byte-unit suffixes accepted by every field in :const:`BYTE_FIELDS`, mapped to the output
#: dtype and the divisor applied to the raw byte count. Multiples are binary (1024-based),
#: matching the convention used elsewhere in the lab module repository.
BYTE_UNITS: dict[str, tuple[pl.DataType, int]] = {
    "b": (pl.Int64, 1),
    "mb": (pl.Float32, 1024 ** 2),
    "gb": (pl.Float32, 1024 ** 3),
}

#: Columns whose name already carries its unit and which are requested verbatim, mapped to
#: their output dtype and help text.
SCALAR_COLUMNS: dict[str, tuple[pl.DataType, str]] = {
    "time": (
        pl.Datetime("us", "UTC"),
        "Wall-clock timestamp of the sample (UTC).",
    ),
    "elapsed_s": (
        pl.Float64,
        "Seconds since the rule started.",
    ),
    "cores_n": (
        pl.Float64,
        "CPU cores used since the previous sample (1.0 is one core fully busy).",
    ),
    "cpu_s": (
        pl.Float64,
        "Cumulative CPU seconds consumed by the tree since the rule started, including reaped "
        "descendants.",
    ),
    "proc_n": (
        pl.Int32,
        "Number of processes in the tree.",
    ),
    "threads_n": (
        pl.Int32,
        "Number of threads across the tree.",
    ),
    "fd_n": (
        pl.Int32,
        "Number of open file descriptors across the tree.",
    ),
    "mem_n": (
        pl.Float64,
        "Tree RSS as a fraction of total system memory (1.0 is all of it).",
    ),
    "read_bps": (
        pl.Float64,
        "Bytes read per second since the previous sample (Linux only; null elsewhere).",
    ),
    "write_bps": (
        pl.Float64,
        "Bytes written per second since the previous sample (Linux only; null elsewhere).",
    ),
}

#: Raw sample fields holding a cumulative counter, reported from the first sample (see :meth:`PerfLog._record`).
CUMULATIVE_FIELDS: tuple[str, ...] = ("cpu_s", "read", "write")

#: Dtype of the constant meta columns. Everything is kept as a string: these are identifiers
#: (rule names, assembly names), not measurements.
META_DTYPE: pl.DataType = pl.Utf8


@dataclass(frozen=True, slots=True)
class Column:
    """A resolved output column: what to name it, where its value comes from, and its type.

    :ivar name: Output column name exactly as it appears in the table (e.g. ``rss_gb``).
    :ivar field: Key of the raw value in a sample dictionary (e.g. ``rss``).
    :ivar dtype: Polars dtype of the output column.
    :ivar divisor: Divisor applied to the raw value; ``1`` leaves it unchanged.
    """

    name: str
    field: str
    dtype: pl.DataType
    divisor: int


def make_column(field: str, unit: Optional[str]) -> Column:
    """Build the output column for a field rendered in a given unit.

    :param field: Sample field name, from :const:`SCALAR_COLUMNS` or :const:`BYTE_FIELDS`.
    :type field: str
    :param unit: Byte unit from :const:`BYTE_UNITS`, or ``None`` for a scalar field.
    :type unit: Optional[str]
    :return: The resolved column.
    :rtype: Column
    """
    if unit is None:
        dtype, _ = SCALAR_COLUMNS[field]
        return Column(name=field, field=field, dtype=dtype, divisor=1)

    dtype, divisor = BYTE_UNITS[unit]
    return Column(name=f"{field}_{unit}", field=field, dtype=dtype, divisor=divisor)


def resolve_field(name: str) -> tuple[str, Optional[str]]:
    """Resolve one requested column name to the field it reads and the unit it is shown in.

    Names in :const:`SCALAR_COLUMNS` are taken verbatim and carry no separate unit. Otherwise
    the name is split on its final underscore into a field from :const:`BYTE_FIELDS` and a unit
    from :const:`BYTE_UNITS`.

    :param name: Requested column name, e.g. ``cores_n`` or ``rss_gb``.
    :type name: str
    :return: Tuple of ``(field, unit)``, where ``unit`` is ``None`` for a scalar column.
    :rtype: tuple[str, Optional[str]]
    :raises ValueError: If the name is not a known column.
    """
    if name in SCALAR_COLUMNS:
        return name, None

    field, _, unit = name.rpartition("_")
    if field in BYTE_FIELDS and unit in BYTE_UNITS:
        return field, unit

    # Point a bare or mis-suffixed byte field at its valid forms before the generic message.
    if name in BYTE_FIELDS or field in BYTE_FIELDS:
        base = name if name in BYTE_FIELDS else field
        forms = ", ".join(f"{base}_{unit}" for unit in BYTE_UNITS)
        raise ValueError(f"Column {name!r} needs a byte unit; use one of: {forms}")

    raise ValueError(f"Unknown column: {name!r}. Available columns: {', '.join(column_names())}")


def resolve_group(name: str) -> tuple[Optional[str], Optional[str]]:
    """Recognize a group keyword and the byte unit it applies to its group.

    ``all`` and ``default`` name a set of columns; a byte-unit suffix on either sets the unit
    every byte column in that set is rendered in. A bare keyword means :const:`DEFAULT_UNIT`.

    :param name: Token from the column list.
    :type name: str
    :return: Tuple of ``(keyword, unit)``, or ``(None, None)`` if this is not a group keyword.
    :rtype: tuple[Optional[str], Optional[str]]
    :raises ValueError: If a group keyword carries a suffix that is not a byte unit.
    """
    if name in GROUP_KEYWORDS:
        return name, DEFAULT_UNIT

    keyword, _, unit = name.rpartition("_")
    if keyword in GROUP_KEYWORDS:
        if unit not in BYTE_UNITS:
            forms = ", ".join(f"{keyword}_{known}" for known in BYTE_UNITS)
            raise ValueError(f"Column {name!r} needs a byte unit; use one of: {keyword}, {forms}")
        return keyword, unit

    return None, None


def column_names() -> list[str]:
    """List every column name and keyword that the column list accepts.

    :return: Group keywords, scalar column names, then every byte field / unit combination.
    :rtype: list[str]
    """
    names = list(GROUP_KEYWORDS)
    names.extend(f"{keyword}_{unit}" for keyword in GROUP_KEYWORDS for unit in BYTE_UNITS)
    names.extend(SCALAR_COLUMNS)
    names.extend(f"{field}_{unit}" for field in BYTE_FIELDS for unit in BYTE_UNITS)
    return names


def parse_columns(spec: str) -> tuple[Column, ...]:
    """Parse and validate a column list.

    The value is a comma-separated list that *adds to* the default columns rather than
    replacing them, so the defaults are always present. A token is either a column name
    (``uss_b``, ``cpu_s``) or one of the group keywords from :const:`GROUP_KEYWORDS`: ``all``
    adds every remaining column, and both keywords accept a byte-unit suffix that sets the unit
    for the byte columns in their group (``all_gb``, ``default_mb``).

    A byte field's unit is taken from the most specific instruction covering it: naming the
    column outright beats ``default_*``, which beats ``all_*``, which beats :const:`DEFAULT_UNIT`.
    Naming one field in several units (``rss_b,rss_gb``) emits a column for each.

    Default columns keep their canonical positions, additions from ``all`` follow in the order
    of :const:`ALL_FIELDS`, and any remaining columns follow in the order requested.

    :param spec: Raw comma-separated column list; empty means the defaults alone.
    :type spec: str
    :return: Resolved columns in output order.
    :rtype: tuple[Column, ...]
    :raises ValueError: If a name is not recognized or a group keyword is given two units.
    """
    group_units: dict[str, str] = {}
    explicit_units: dict[str, list[Optional[str]]] = {}
    requested: list[str] = []

    for token in spec.split(","):
        token = token.strip().lower()
        if not token:
            continue

        keyword, unit = resolve_group(token)
        if keyword is not None:
            previous = group_units.get(keyword)
            if previous is not None and previous != unit:
                raise ValueError(
                    f"Conflicting units for {keyword!r}: {keyword}_{previous} and {keyword}_{unit}"
                )
            group_units[keyword] = unit
            continue

        field, unit = resolve_field(token)
        units = explicit_units.setdefault(field, [])
        if unit not in units:
            units.append(unit)
        if field not in requested:
            requested.append(field)

    # The defaults always lead, then everything "all" pulls in, then anything named directly.
    fields: list[str] = list(DEFAULT_FIELDS)
    if "all" in group_units:
        fields.extend(field for field in ALL_FIELDS if field not in fields)
    fields.extend(field for field in requested if field not in fields)

    columns: list[Column] = []

    for field in fields:
        if field in SCALAR_COLUMNS:
            columns.append(make_column(field, None))
            continue

        units = explicit_units.get(field)
        if units:
            # Named outright, in every unit it was named in.
            columns.extend(make_column(field, unit) for unit in units)
            continue

        # Otherwise the most specific group keyword covering this field decides its unit.
        unit = group_units.get("all", DEFAULT_UNIT)
        if field in DEFAULT_FIELDS:
            unit = group_units.get("default", unit)
        columns.append(make_column(field, unit))

    return tuple(columns)


def apply_meta(frame: pl.DataFrame, meta: dict[str, str]) -> pl.DataFrame:
    """Prepend the constant meta columns to the sample table.

    :param frame: Table of samples.
    :type frame: pl.DataFrame
    :param meta: Mapping of column name to constant value.
    :type meta: dict[str, str]
    :return: The table with the meta columns first, unchanged if there are none.
    :rtype: pl.DataFrame
    """
    if not meta:
        return frame

    # Leading columns, so a row identifies itself before its measurements. A literal
    # broadcasts over however many rows there are, including none.
    return frame.with_columns(
        pl.lit(value, dtype=META_DTYPE).alias(key) for key, value in meta.items()
    ).select(list(meta) + frame.columns)


@dataclass(slots=True)
class Accumulator:
    """Per-sample running totals gathered while walking the process tree.

    :ivar proc_n: Processes visited.
    :ivar threads_n: Threads summed over the tree.
    :ivar fd_n: Open file descriptors summed over the tree, ``None`` if never readable.
    :ivar rss: Resident bytes summed over the tree.
    :ivar vss: Virtual bytes summed over the tree.
    :ivar pss: Proportional-set-size bytes, ``None`` if unsupported on this platform.
    :ivar uss: Unique-set-size bytes, ``None`` if unsupported on this platform.
    :ivar swap: Swapped-out bytes, ``None`` if unsupported on this platform.
    :ivar cpu_own: CPU seconds consumed by the live processes themselves.
    :ivar cpu_children: CPU seconds of descendants the live processes have already reaped.
    :ivar cpu_by_key: Own CPU seconds of each live process, keyed by ``(pid, create_time)``.
    :ivar io_by_key: ``(read_bytes, write_bytes)`` of each live process, keyed the same way.
    :ivar io_seen: Whether any process reported I/O counters.
    """

    proc_n: int = 0
    threads_n: int = 0
    fd_n: Optional[int] = None
    rss: int = 0
    vss: int = 0
    pss: Optional[int] = None
    uss: Optional[int] = None
    swap: Optional[int] = None
    cpu_own: float = 0.0
    cpu_children: float = 0.0
    cpu_by_key: dict[tuple[int, float], float] = field(default_factory=dict)
    io_by_key: dict[tuple[int, float], tuple[int, int]] = field(default_factory=dict)
    io_seen: bool = False


class TreeMonitor:
    """Samples the resource usage of a process and every descendant it still owns.

    One instance follows one process tree for the life of the monitored rule. It keeps the
    running counters needed to turn per-process cumulative values (CPU seconds, I/O bytes) into
    per-interval rates, and holds those counters monotonic so that a process disappearing
    between two samples cannot produce a negative rate.

    :ivar _root: The monitored process.
    :ivar _start: ``time.monotonic()`` reading taken when monitoring started.
    :ivar _fields: Set of raw sample fields the selected columns actually need; expensive
        readings are skipped when nothing asks for them.
    :ivar _total_memory: Total system memory in bytes, used for ``mem_n``.
    :ivar _cpu_s: Monotonic cumulative CPU seconds consumed by the tree.
    :ivar _read_b: Monotonic cumulative bytes read by the tree, ``None`` if unsupported.
    :ivar _write_b: Monotonic cumulative bytes written by the tree, ``None`` if unsupported.
    :ivar _live_cpu: Own CPU seconds of the processes alive at the previous sample.
    :ivar _retired_cpu_s: Own CPU seconds of processes seen alive that have since exited.
    :ivar _live_io: I/O counters of the processes alive at the previous sample.
    :ivar _retired_read_b: Bytes read by processes seen alive that have since exited.
    :ivar _retired_write_b: Bytes written by processes seen alive that have since exited.
    :ivar _prev_monotonic: ``time.monotonic()`` reading of the previous sample.
    :ivar _prev_cpu_s: Cumulative CPU seconds at the previous sample.
    :ivar _prev_read_b: Cumulative bytes read at the previous sample.
    :ivar _prev_write_b: Cumulative bytes written at the previous sample.
    :ivar _first: Whether the next sample is the first one (no interval to rate against).
    :ivar _warned: Readings that have already produced a one-time warning.
    """

    def __init__(self, pid: int, fields: frozenset[str], start_monotonic: float) -> None:
        """Start following the tree rooted at ``pid``.

        :param pid: Process ID of the root of the tree.
        :type pid: int
        :param fields: Raw sample fields the selected columns need.
        :type fields: frozenset[str]
        :param start_monotonic: ``time.monotonic()`` taken when monitoring started.
        :type start_monotonic: float
        """
        self._root: psutil.Process = psutil.Process(pid)
        self._start: float = start_monotonic
        self._fields: frozenset[str] = fields
        self._total_memory: int = psutil.virtual_memory().total

        self._cpu_s: float = 0.0
        self._read_b: Optional[int] = None
        self._write_b: Optional[int] = None

        # Values carried by processes that were alive at the previous sample, so that what a
        # process consumed before exiting can be retired into the cumulative totals rather
        # than vanishing. Keyed by (pid, create_time) so a recycled PID is never mistaken for
        # the process it replaced.
        self._live_cpu: dict[tuple[int, float], float] = {}
        self._retired_cpu_s: float = 0.0
        self._live_io: dict[tuple[int, float], tuple[int, int]] = {}
        self._retired_read_b: int = 0
        self._retired_write_b: int = 0

        self._prev_monotonic: float = start_monotonic
        self._prev_cpu_s: float = 0.0
        self._prev_read_b: int = 0
        self._prev_write_b: int = 0

        self._first: bool = True
        self._warned: set[str] = set()

    def _processes(self) -> list[psutil.Process]:
        """Collect the root process and all of its live descendants.

        :return: The tree's processes, or an empty list once the root is gone.
        :rtype: list[psutil.Process]
        """
        try:
            return [self._root] + self._root.children(recursive=True)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            return []

    def _warn_once(self, key: str, message: str, *args: Any) -> None:
        """Emit a warning the first time a given reading fails.

        Readings can fail for the whole run (a platform that does not implement them) or for
        individual processes, and either way the warning is worth exactly one line.

        :param key: Identifier of the reading, used to suppress repeats.
        :type key: str
        :param message: ``%``-style log message.
        :type message: str
        :param args: Arguments interpolated into the message.
        :type args: Any
        """
        if key not in self._warned:
            self._warned.add(key)
            logger.warning(message, *args)

    def _read_process(self, proc: psutil.Process, acc: Accumulator) -> None:
        """Read one process into the running totals for this sample.

        Memory and CPU are required; every other reading is optional and independently guarded
        so that one unavailable field does not discard the rest of the process's contribution.

        :param proc: The process to read.
        :type proc: psutil.Process
        :param acc: Running totals for this sample, updated in place.
        :type acc: Accumulator
        :raises psutil.NoSuchProcess: If the process exits while being read.
        :raises psutil.AccessDenied: If its basic memory or CPU counters cannot be read.
        """
        key = (proc.pid, proc.create_time())

        memory = proc.memory_info()
        times = proc.cpu_times()

        acc.proc_n += 1
        acc.rss += memory.rss
        acc.vss += memory.vms

        own = times.user + times.system
        acc.cpu_own += own
        acc.cpu_by_key[key] = own

        # CPU of descendants this process has already reaped. Reported on Linux, 0.0 elsewhere.
        acc.cpu_children += getattr(times, "children_user", 0.0)
        acc.cpu_children += getattr(times, "children_system", 0.0)

        if "threads_n" in self._fields:
            try:
                acc.threads_n += proc.num_threads()
            except psutil.AccessDenied:
                self._warn_once("threads", "Thread counts unavailable for some processes")

        if "fd_n" in self._fields:
            try:
                acc.fd_n = (acc.fd_n or 0) + proc.num_fds()
            except (psutil.AccessDenied, AttributeError, NotImplementedError):
                self._warn_once("fd", "File descriptor counts unavailable; fd_n may be low")

        if self._fields & {"pss", "uss", "swap"}:
            try:
                full = proc.memory_full_info()
                # pss and swap come from Linux smaps; uss is available more widely.
                if (value := getattr(full, "pss", None)) is not None:
                    acc.pss = (acc.pss or 0) + value
                if (value := getattr(full, "uss", None)) is not None:
                    acc.uss = (acc.uss or 0) + value
                if (value := getattr(full, "swap", None)) is not None:
                    acc.swap = (acc.swap or 0) + value
            except (psutil.AccessDenied, NotImplementedError):
                self._warn_once(
                    "full",
                    "Detailed memory (pss/uss/swap) unavailable on this platform or for these "
                    "processes; those columns will be null. rss/vss are unaffected.",
                )

        if self._fields & {"read", "write", "read_bps", "write_bps"}:
            try:
                counters = proc.io_counters()
                acc.io_by_key[key] = (counters.read_bytes, counters.write_bytes)
                acc.io_seen = True
            except (psutil.AccessDenied, AttributeError, NotImplementedError):
                self._warn_once(
                    "io",
                    "I/O counters unavailable on this platform or for these processes; the "
                    "read/write columns will be null.",
                )

    def sample(self) -> dict[str, Any]:
        """Take one snapshot of the tree.

        :return: Raw sample values keyed by field name. Fields whose source is unsupported on
            this platform, or which no selected column needs, are ``None``.
        :rtype: dict[str, Any]
        """
        now_monotonic = time.monotonic()
        timestamp = datetime.now(timezone.utc)
        interval_s = now_monotonic - self._prev_monotonic

        acc = Accumulator()

        for proc in self._processes():
            try:
                # oneshot() batches the underlying per-process reads into one system call
                # where the platform allows it.
                with proc.oneshot():
                    self._read_process(proc, acc)
            except psutil.NoSuchProcess:
                # The process exited while being read; it drops out of this sample.
                continue
            except psutil.AccessDenied:
                self._warn_once(
                    "core",
                    "Some processes could not be inspected at all; reported usage may be low",
                )
                continue

        self._retire(acc)

        # Two independent lower bounds on the tree's cumulative CPU, neither of which double
        # counts. "Reaped" covers descendants that lived and died entirely between two samples
        # but is 0.0 where children_* is not reported; "observed" covers every process a sample
        # actually saw. The larger is the better estimate, and max() with the previous value
        # keeps the counter monotonic while a process sits unreaped.
        by_reaped = acc.cpu_own + acc.cpu_children
        by_observed = acc.cpu_own + self._retired_cpu_s
        self._cpu_s = max(self._cpu_s, by_reaped, by_observed)

        if acc.io_seen:
            read_total = sum(read for read, _ in acc.io_by_key.values()) + self._retired_read_b
            write_total = sum(write for _, write in acc.io_by_key.values()) + self._retired_write_b
            self._read_b = max(self._read_b or 0, read_total)
            self._write_b = max(self._write_b or 0, write_total)

        # The first sample has no preceding interval to difference against, so its rates are
        # reported as zero rather than as a cumulative total divided by the process start-up
        # time, which would read as a meaningless spike.
        if self._first:
            cores_n = 0.0
            read_bps = None if self._read_b is None else 0.0
            write_bps = None if self._write_b is None else 0.0
        else:
            cores_n = self._rate(self._cpu_s, self._prev_cpu_s, interval_s)
            read_bps = (
                None if self._read_b is None
                else self._rate(self._read_b, self._prev_read_b, interval_s)
            )
            write_bps = (
                None if self._write_b is None
                else self._rate(self._write_b, self._prev_write_b, interval_s)
            )

        self._first = False
        self._prev_monotonic = now_monotonic
        self._prev_cpu_s = self._cpu_s
        self._prev_read_b = self._read_b or 0
        self._prev_write_b = self._write_b or 0

        return {
            "time": timestamp,
            "elapsed_s": now_monotonic - self._start,
            "pss": acc.pss,
            "rss": acc.rss,
            "vss": acc.vss,
            "uss": acc.uss,
            "swap": acc.swap,
            "read": self._read_b,
            "write": self._write_b,
            "cores_n": cores_n,
            "cpu_s": self._cpu_s,
            "proc_n": acc.proc_n,
            "threads_n": acc.threads_n,
            "fd_n": acc.fd_n,
            "mem_n": acc.rss / self._total_memory if self._total_memory else 0.0,
            "read_bps": read_bps,
            "write_bps": write_bps,
        }

    def _retire(self, acc: Accumulator) -> None:
        """Fold processes that have exited since the previous sample into the retired totals.

        Only each process's *own* consumption is retired. Retiring the ``children_*`` figures
        as well would double-count any descendant that a sample had already observed alive.

        :param acc: Running totals for this sample, holding the currently live processes.
        :type acc: Accumulator
        """
        for key, value in self._live_cpu.items():
            if key not in acc.cpu_by_key:
                self._retired_cpu_s += value
        self._live_cpu = acc.cpu_by_key

        for key, (read, write) in self._live_io.items():
            if key not in acc.io_by_key:
                self._retired_read_b += read
                self._retired_write_b += write
        self._live_io = acc.io_by_key

    @staticmethod
    def _rate(current: float, previous: float, interval_s: float) -> float:
        """Convert the change in a cumulative counter into a per-second rate.

        :param current: Counter value at this sample.
        :type current: float
        :param previous: Counter value at the previous sample.
        :type previous: float
        :param interval_s: Wall-clock seconds between the two samples.
        :type interval_s: float
        :return: Change per second, or ``0.0`` if no time has passed.
        :rtype: float
        """
        if interval_s <= 0:
            return 0.0
        return (current - previous) / interval_s


def build_schema(columns: tuple[Column, ...]) -> dict[str, pl.DataType]:
    """Build the output schema for the selected columns.

    :param columns: Resolved columns in output order.
    :type columns: tuple[Column, ...]
    :return: Mapping of column name to dtype, in output order.
    :rtype: dict[str, pl.DataType]
    """
    return {column.name: column.dtype for column in columns}


def build_row(sample: dict[str, Any], columns: tuple[Column, ...]) -> dict[str, Any]:
    """Project one raw sample onto the selected columns, applying unit conversions.

    :param sample: Raw sample from :meth:`TreeMonitor.sample`.
    :type sample: dict[str, Any]
    :param columns: Resolved columns in output order.
    :type columns: tuple[Column, ...]
    :return: One output row keyed by output column name.
    :rtype: dict[str, Any]
    """
    row: dict[str, Any] = {}

    for column in columns:
        value = sample[column.field]
        if value is not None and column.divisor != 1:
            value = value / column.divisor
        row[column.name] = value

    return row


class PerfLog:
    """Sample the resource usage of this process and its descendants while a block of code runs.

    The first sample is taken on entry, then one every ``interval_s`` seconds on a fixed schedule so sampling does not
    drift; if a sample is delayed past the next scheduled time, the missed slots are skipped rather than fired
    back-to-back. A final sample is taken on exit. The table is written to ``path`` after every sample and on exit,
    including when the block raises an exception.

    Use as a context manager::

        with PerfLog('log/perf/HG00171/call_inter_h1.parquet', 30, meta={'rule': 'call_inter'}):
            ...

    :ivar _path: Output parquet file.
    :ivar _interval_s: Sampling interval in seconds.
    :ivar _columns: Resolved output columns.
    :ivar _schema: Output schema of the sampled columns.
    :ivar _meta: Constant columns written ahead of the sampled ones.
    :ivar _pid: Process ID at the root of the sampled tree.
    :ivar _rows: Samples collected so far.
    :ivar _baseline: Cumulative counters at the first sample, subtracted from every sample.
    :ivar _monitor: Process tree monitor, set on entry.
    :ivar _start: ``time.monotonic()`` reading taken on entry.
    :ivar _stop: Set on exit to stop the sampling thread.
    :ivar _thread: Sampling thread, set on entry.
    :ivar _failed: Set if sampling failed and stopped early.
    :ivar _write_warned: Set once a failed write has been reported.
    """

    def __init__(
            self,
            path: Path | str,
            interval_s: float,
            columns: str = PERF_COLUMNS,
            meta: Optional[dict[str, str]] = None,
            pid: Optional[int] = None,
    ) -> None:
        """Prepare to sample.

        :param path: Output parquet file. Parent directories are created on entry.
        :param interval_s: Sampling interval in seconds.
        :param columns: Comma-separated columns to add to the default columns (see :func:`parse_columns`).
        :param meta: Constant columns written ahead of the sampled ones in every row.
        :param pid: Process ID at the root of the sampled tree. Defaults to this process.

        :raises ValueError: If the interval is not positive, the column list is not valid, or a meta column collides
            with a sampled column.
        """
        if interval_s <= 0:
            raise ValueError(f'Performance sampling interval must be positive: {interval_s}')

        self._path: Path = Path(path)
        self._interval_s: float = float(interval_s)
        self._columns: tuple[Column, ...] = parse_columns(columns)
        self._schema: dict[str, pl.DataType] = build_schema(self._columns)
        self._meta: dict[str, str] = {key: str(value) for key, value in (meta or {}).items()}
        self._pid: int = os.getpid() if pid is None else pid

        # A meta column shadowing a sampled one would silently discard the measurement.
        for key in self._meta:
            if key in self._schema or key in column_names():
                raise ValueError(f'Performance log meta column {key!r} collides with a sampled column')

        self._rows: list[dict[str, Any]] = []
        self._baseline: Optional[dict[str, Any]] = None
        self._monitor: Optional[TreeMonitor] = None
        self._start: Optional[float] = None
        self._stop: threading.Event = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self._failed: bool = False
        self._write_warned: bool = False

    @property
    def path(self) -> Path:
        """Output parquet file."""
        return self._path

    def table(self) -> pl.DataFrame:
        """Get the samples collected so far.

        :returns: Sample table with meta columns first.
        """
        return apply_meta(pl.DataFrame(self._rows, schema=self._schema), self._meta)

    def __enter__(self) -> 'PerfLog':
        """Take the first sample and start the sampling thread.

        :returns: This object.
        """
        self._path.parent.mkdir(parents=True, exist_ok=True)

        self._start = time.monotonic()
        self._monitor = TreeMonitor(
            self._pid, frozenset(column.field for column in self._columns), self._start
        )

        # The first sample is written before the rule does any work, so a table that cannot be written fails the rule
        # at the start instead of after it has run for hours.
        self._record()
        self._write()

        self._thread = threading.Thread(target=self._run, name=f'perf_log:{self._path.name}', daemon=True)
        self._thread.start()

        return self

    def __exit__(
            self,
            exc_type: Optional[type[BaseException]],
            exc_value: Optional[BaseException],
            traceback: Optional[TracebackType],
    ) -> None:
        """Stop sampling, take a final sample, and write the table.

        :param exc_type: Exception type raised by the block, if any.
        :param exc_value: Exception raised by the block, if any.
        :param traceback: Traceback of the exception, if any.

        :raises OSError: If the table cannot be written and the block did not raise an exception. If it did, the
            failure is logged so the block's own exception is the one raised.
        """
        self._stop.set()

        if self._thread is not None:
            self._thread.join()

        if not self._failed:
            try:
                self._record()
            except Exception:
                logger.warning('Failed to take the final performance sample: %s', self._path, exc_info=True)

        try:
            self._write()
        except Exception:
            if exc_type is None:
                raise

            logger.warning('Failed to write performance log: %s', self._path, exc_info=True)

    def _record(self) -> None:
        """Take one sample and add it to the table.

        Cumulative counters are reported from the start of the rule. The tree is rooted at the Snakemake job process,
        which spent CPU and I/O starting Snakemake before the rule started, so the first sample is the baseline.
        """
        sample = self._monitor.sample()

        if self._baseline is None:
            self._baseline = {key: sample[key] for key in CUMULATIVE_FIELDS}

        for key, base in self._baseline.items():
            if sample[key] is not None and base is not None:
                sample[key] -= base

        self._rows.append(build_row(sample, self._columns))

    def _write(self) -> None:
        """Write the table, replacing the previous one atomically.

        A job killed mid-write keeps the previous table instead of a truncated file. The temporary file is hidden so a
        glob for ``*.parquet`` never picks it up.
        """
        temp_path = self._path.with_name(f'.{self._path.name}.tmp')

        self.table().write_parquet(temp_path)
        os.replace(temp_path, self._path)

    def _run(self) -> None:
        """Sample on a fixed schedule until stopped."""
        tick = 1

        try:
            while True:
                delay = self._start + tick * self._interval_s - time.monotonic()

                if delay > 0 and self._stop.wait(delay):
                    break

                self._record()

                # A failed write leaves the previous table in place, and the next write includes every sample, so a
                # transient failure loses nothing. The write on exit reports a failure that persists.
                try:
                    self._write()
                except OSError:
                    if not self._write_warned:
                        self._write_warned = True
                        logger.warning('Failed to write performance log: %s', self._path, exc_info=True)

                # Advance past any slots missed while sampling, keeping the schedule absolute.
                now = time.monotonic()
                tick += 1
                while self._start + tick * self._interval_s <= now:
                    tick += 1

        except Exception:
            self._failed = True
            logger.warning(
                'Performance sampling stopped early; the samples collected so far are kept: %s', self._path,
                exc_info=True
            )


def perf_log_path(
        rule: str,
        wildcards: Any,
) -> Path:
    """Get the performance log path for a rule.

    The path is "log/perf/{asm_name}/{rule}_{wildcards}.parquet" where "wildcards" is the value of each wildcard
    except "asm_name" joined by "_" in the order they appear in the rule. Rules without an "asm_name" wildcard write to
    "log/perf/" directly.

    :param rule: Rule name.
    :param wildcards: Rule wildcards (Snakemake wildcards object or a dict).

    :returns: Path to the performance log, relative to the PAV run directory.
    """
    wildcards = dict(wildcards.items())

    asm_name = wildcards.pop('asm_name', None)

    file_name = '_'.join([rule] + [str(val) for val in wildcards.values()]) + '.parquet'

    return (PERF_LOG_DIR / asm_name if asm_name is not None else PERF_LOG_DIR) / file_name


def perf_log(
        rule: str,
        wildcards: Any,
        pav_config: Optional[dict] = None,
        asm_table: Optional[pl.DataFrame] = None,
) -> contextlib.AbstractContextManager:
    """Get a context manager logging the performance of a rule if config parameter "perf_log" enables it.

    In a rule::

        run:
            with pav3.workflow.perf.perf_log(rule, wildcards, PAV_CONFIG, ASM_TABLE):
                ...

    :param rule: Rule name.
    :param wildcards: Rule wildcards (Snakemake wildcards object or a dict).
    :param pav_config: PAV configuration.
    :param asm_table: Assembly table. The assembly's CONFIG column may override "perf_log".

    :returns: A :class:`PerfLog` writing to :func:`perf_log_path`, or a context manager doing nothing if performance
        logging is disabled.
    """
    wildcards = dict(wildcards.items())

    asm_name = wildcards.get('asm_name', None)

    interval_s = PavParams(asm_name, pav_config, asm_table).perf_log

    if interval_s == 0:
        return contextlib.nullcontext()

    return PerfLog(
        perf_log_path(rule, wildcards),
        interval_s,
        meta={
            'rule': rule,
            'asm_name': asm_name if asm_name is not None else '',
            'hap': wildcards.get('hap', ''),
            'wildcards': ','.join(f'{key}={val}' for key, val in wildcards.items() if key != 'asm_name'),
        }
    )
