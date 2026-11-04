"""Shared table and FASTA input/output for the probe tools.

CSV is the default table format, ``.tsv``/``.txt`` use tabs, and ``.xlsx``
uses Excel. Identifier and sequence cells are read as text, including values
such as ``001`` and ``NA``. Empty cells stay empty instead of becoming ``nan``.
"""

import os
from pathlib import Path
import tempfile

import pandas as pd


def read_probe_table(path, required_columns=()):
    """Read a table without changing identifiers or inventing missing sequences.

    Other columns retain pandas' normal numeric/boolean inference. Callers may
    require specific columns; the reverse-complement tool can read any schema.
    A table with a header and no data rows is allowed.
    """
    source = Path(path)
    options = {"dtype": {"id": str, "sequence": str}, "keep_default_na": False}
    try:
        if source.suffix.lower() == ".xlsx":
            frame = pd.read_excel(source, engine="openpyxl", **options)
        else:
            separator = "\t" if source.suffix.lower() in {".tsv", ".txt"} else ","
            frame = pd.read_csv(source, sep=separator, encoding="utf-8-sig", **options)
    except pd.errors.EmptyDataError as exc:
        raise ValueError(f"输入文件为空或没有表头: {source}") from exc

    if len(frame.columns) == 0:
        raise ValueError(f"输入文件为空或没有表头: {source}")
    if isinstance(required_columns, str):
        required_columns = (required_columns,)
    missing = [column for column in required_columns if column not in frame.columns]
    if missing:
        raise ValueError(f"输入文件缺少必要的列: {', '.join(missing)}")

    for column in ("id", "sequence"):
        if column in frame.columns:
            frame[column] = frame[column].map(
                lambda value: "" if pd.isna(value) else str(value)
            )
    return frame


def write_probe_table(frame, path):
    """Write a table atomically in the format selected by its extension.

    A failed write leaves any previous destination untouched. Temporary files
    are created beside the destination so replacement stays on one filesystem.
    """
    destination = Path(path)
    if not destination.name:
        raise ValueError("请选择输出文件")
    if len(frame.columns) == 0:
        raise ValueError("输出数据没有列名")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        prefix=f".{destination.name}.",
        suffix=destination.suffix or ".csv",
        dir=destination.parent,
        delete=False,
    ) as temporary:
        temporary_path = Path(temporary.name)
    try:
        if destination.suffix.lower() == ".xlsx":
            frame.to_excel(temporary_path, index=False, engine="openpyxl")
        else:
            separator = "\t" if destination.suffix.lower() in {".tsv", ".txt"} else ","
            frame.to_csv(temporary_path, sep=separator, index=False, encoding="utf-8-sig")
        os.replace(temporary_path, destination)
    finally:
        temporary_path.unlink(missing_ok=True)


def read_fasta_records(path):
    """Return separate ``(identifier, sequence)`` records from a FASTA file.

    Sequence whitespace is removed and letters are uppercased. Non-whitespace
    characters are retained for the caller's alphabet validation, never silently
    discarded. Empty records and sequence text before the first header are
    rejected. Repeated record identifiers remain separate records.
    """
    records = []
    identifier = None
    fragments = []

    def finish_record():
        sequence = "".join(fragments).upper()
        if not sequence:
            raise ValueError(f"FASTA记录 {identifier} 没有序列")
        records.append((identifier, sequence))

    with open(path, "r", encoding="utf-8-sig") as source:
        for line_number, line in enumerate(source, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    finish_record()
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"FASTA第 {line_number} 行缺少记录名称")
                identifier = header.split()[0]
                fragments = []
            else:
                if identifier is None:
                    raise ValueError(f"FASTA第 {line_number} 行的序列前缺少 >记录名称")
                fragments.append("".join(line.split()))

    if identifier is None:
        raise ValueError("FASTA文件中没有序列记录")
    finish_record()
    return records
