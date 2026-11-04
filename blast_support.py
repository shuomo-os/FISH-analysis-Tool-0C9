"""Shared BLAST operations for the analysis and design tabs.

The reported counts describe the selected database and search settings, not a
biological specificity score. Multiple subject accessions need not be different
genes. See README.md for the original project's attribution and license.

Short-query settings follow NCBI's BLAST Program Selection Guide, section 4.3:
https://blast.ncbi.nlm.nih.gov/blast/BLAST_guide.pdf
The tabular field definitions are documented in the BLAST+ manual:
https://www.ncbi.nlm.nih.gov/books/NBK279684/
"""

import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from numbers import Integral


MIN_QUERY_LENGTH = 7
DEFAULT_MAX_TARGET_SEQS = 1000
DEFAULT_EVALUE = 1000.0
BLAST_FIELDS = (
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "qcovhsp",
)
BLAST_OUTFMT = "7 " + " ".join(BLAST_FIELDS)


class BlastError(RuntimeError):
    """A BLAST command or its output could not be used."""


class BlastCancelled(BlastError):
    """The user cancelled the running search."""


def write_query_fasta(queries, path):
    """Write ``(source_row_position, DNA_sequence)`` pairs using safe query IDs.

    Return an internal-ID-to-source-row-position dictionary. The user's own ID
    is intentionally not used: it may contain spaces, leading zeros, or
    duplicates. Validate all input before opening the destination.
    """
    records = []
    query_rows = {}
    seen_positions = set()
    for number, (position, sequence) in enumerate(queries, 1):
        if (not isinstance(position, Integral) or isinstance(position, bool)
                or position < 0 or position in seen_positions):
            raise ValueError("BLAST 查询的原始行位置必须是互不重复的非负整数。")
        if not isinstance(sequence, str):
            raise ValueError(f"第 {position + 1} 行的 BLAST 序列必须是文本。")
        sequence = sequence.strip().upper()
        if not sequence or set(sequence) - set("ACGT"):
            raise ValueError(f"第 {position + 1} 行的 BLAST 序列只允许 A、C、G、T。")
        if len(sequence) < MIN_QUERY_LENGTH:
            raise ValueError(
                f"第 {position + 1} 行的序列短于 {MIN_QUERY_LENGTH} nt，"
                "不能使用当前 BLAST 种子长度进行检索。"
            )
        query_id = f"q{number:06d}"
        records.append(f">{query_id}\n{sequence}\n")
        query_rows[query_id] = int(position)
        seen_positions.add(position)
    if not records:
        raise ValueError("没有可供 BLAST 检索的序列。")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(records)
    return query_rows


def build_blast_command(executable, database, query_path, output_path,
                        max_target_seqs=DEFAULT_MAX_TARGET_SEQS,
                        evalue=DEFAULT_EVALUE):
    """Build one explicit short-oligonucleotide search command for both tabs.

    ``max_target_seqs`` limits the number of reported subject sequences. Counts
    derived from this output are not an exhaustive count of database matches.
    """
    if (not isinstance(max_target_seqs, Integral)
            or isinstance(max_target_seqs, bool) or max_target_seqs < 1):
        raise ValueError("BLAST 报告的最大目标序列数必须是正整数。")
    try:
        evalue = float(evalue)
    except (TypeError, ValueError) as exc:
        raise ValueError("BLAST E-value 必须是正的有限数值。") from exc
    if not math.isfinite(evalue) or evalue <= 0:
        raise ValueError("BLAST E-value 必须是正的有限数值。")
    paths = [os.fspath(value) for value in
             (executable, database, query_path, output_path)]
    if any(not isinstance(value, str) or not value.strip() for value in paths):
        raise ValueError("BLAST 程序、数据库、查询和输出路径不能为空。")
    executable, database, query_path, output_path = paths
    # NCBI parses -db as a database list even after argv handling. Quote this
    # application's single database prefix internally when its path has spaces.
    if any(character.isspace() for character in database):
        if not (database.startswith('"') and database.endswith('"')):
            if '"' in database:
                raise ValueError("BLAST 数据库前缀不能包含未配对的双引号。")
            database = '"' + database + '"'
    return [
        executable, "-task", "blastn-short", "-word_size", str(MIN_QUERY_LENGTH),
        "-db", database, "-query", query_path, "-out", output_path,
        "-outfmt", BLAST_OUTFMT, "-max_target_seqs", str(max_target_seqs),
        "-evalue", format(evalue, ".15g"), "-dust", "no",
        "-soft_masking", "false", "-strand", "both", "-parse_deflines",
    ]


def _stop_process(process, termination_timeout):
    """Terminate, reap, and drain the child; escalate if it ignores termination."""
    if process.poll() is None:
        try:
            process.terminate()
        except OSError:
            # On Windows a process that just exited may raise a generic OSError.
            if process.poll() is None:
                raise
    try:
        return process.communicate(timeout=termination_timeout)
    except subprocess.TimeoutExpired:
        try:
            process.kill()
        except OSError:
            if process.poll() is None:
                raise
        return process.communicate()


def run_blast(command, cancel_event=None, *, poll_interval=0.2,
              termination_timeout=2.0):
    """Run an argument list without a shell, responding to a cancellation Event.

    Run this function in a worker thread. It deliberately never accesses Tk.
    Return ``subprocess.CompletedProcess`` on success; raise ``BlastError`` on
    failure or ``BlastCancelled`` after the cancelled child has been reaped.
    """
    if isinstance(command, (str, bytes)) or not command:
        raise ValueError("BLAST 命令必须是非空参数列表。")
    if (not math.isfinite(poll_interval) or poll_interval <= 0
            or not math.isfinite(termination_timeout) or termination_timeout <= 0):
        raise ValueError("进程检查间隔和终止超时必须是正的有限数值。")
    if cancel_event is not None and cancel_event.is_set():
        raise BlastCancelled("BLAST 已取消。")
    try:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, encoding="utf-8", errors="replace",
            creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0,
        )
    except OSError as exc:
        raise BlastError(f"无法启动 BLAST：{exc}") from exc
    try:
        while True:
            if cancel_event is not None and cancel_event.is_set():
                _stop_process(process, termination_timeout)
                raise BlastCancelled("BLAST 已取消，检索进程已停止。")
            try:
                stdout, stderr = process.communicate(timeout=poll_interval)
                break
            except subprocess.TimeoutExpired:
                continue
        if cancel_event is not None and cancel_event.is_set():
            raise BlastCancelled("BLAST 已取消。")
        if process.returncode != 0:
            detail = stderr.strip() or stdout.strip() or "程序未返回错误详情。"
            raise BlastError(f"BLAST 退出码 {process.returncode}：{detail}")
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
    finally:
        if process.poll() is None:
            _stop_process(process, termination_timeout)


def publish_blast_report(raw_source, mapping_source, output_path):
    """Publish a validated report and its ID mapping, rolling back on failure.

    Callers must finish parsing the staged report and prevent concurrent writers
    to these paths first. Each individual replacement is atomic. If publishing
    the raw report fails, restore the previous mapping so it cannot be paired
    with the wrong report. Retain its backup if rollback itself fails.
    """
    output = Path(output_path)
    mapping_output = Path(str(output) + ".queries.csv")
    backup = None
    mapping_published = False
    retain_backup = False
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        if os.path.lexists(mapping_output):
            descriptor, backup_name = tempfile.mkstemp(
                prefix=f".{output.name}.queries.backup-", suffix=".csv", dir=output.parent
            )
            os.close(descriptor)
            backup = Path(backup_name)
            # Copy rather than move: a backup failure must leave both outputs intact.
            shutil.copy2(mapping_output, backup)
        os.replace(mapping_source, mapping_output)
        mapping_published = True
        os.replace(raw_source, output)
    except OSError as publication_error:
        if mapping_published:
            try:
                if backup is not None:
                    os.replace(backup, mapping_output)
                else:
                    mapping_output.unlink(missing_ok=True)
            except OSError as rollback_error:
                retain_backup = True
                recovery = (f"旧映射备份保留在: {backup}" if backup is not None else
                            f"该映射原先不存在，请移除本次生成的映射: {mapping_output}")
                raise BlastError(
                    f"BLAST报告发布失败: {publication_error}；原始BLAST报告未替换。"
                    f"查询映射回滚失败: {rollback_error}。{recovery}。"
                    "当前报告与映射可能不对应，请勿合并使用。"
                ) from publication_error
        raise BlastError(
            f"BLAST报告发布失败: {publication_error}；此前的报告和查询映射已保留。"
        ) from publication_error
    finally:
        if backup is not None and not retain_backup:
            try:
                backup.unlink(missing_ok=True)
            except OSError:
                # A leftover backup does not invalidate a correctly published pair.
                pass


def _query_id(value, expected_ids, line_number):
    """Accept the ordinary local ID and NCBI's explicit local-ID spelling."""
    if value in expected_ids:
        return value
    if value.startswith("lcl|") and value[4:] in expected_ids:
        return value[4:]
    raise BlastError(f"BLAST 第 {line_number} 行出现未知查询 ID：{value}")


def _parse_hit(parts, line_number):
    if len(parts) != len(BLAST_FIELDS):
        raise BlastError(
            f"BLAST 第 {line_number} 行应有 {len(BLAST_FIELDS)} 个制表符分隔字段，"
            f"实际为 {len(parts)}；请使用本程序生成的 BLAST 输出。"
        )
    try:
        hit = {
            "subject": parts[1],
            "identity": float(parts[2]),
            "alignment_length": int(parts[3]),
            "mismatches": int(parts[4]),
            "gap_opens": int(parts[5]),
            "query_start": int(parts[6]),
            "query_end": int(parts[7]),
            "subject_start": int(parts[8]),
            "subject_end": int(parts[9]),
            "evalue": float(parts[10]),
            "bitscore": float(parts[11]),
            "query_length": int(parts[12]),
            "query_coverage": float(parts[13]),
        }
    except ValueError as exc:
        raise BlastError(f"BLAST 第 {line_number} 行含无效数值。") from exc
    floats = (hit["identity"], hit["evalue"], hit["bitscore"], hit["query_coverage"])
    if not all(math.isfinite(value) for value in floats):
        raise BlastError(f"BLAST 第 {line_number} 行含非有限数值。")
    query_span = abs(hit["query_end"] - hit["query_start"]) + 1
    subject_span = abs(hit["subject_end"] - hit["subject_start"]) + 1
    if (not hit["subject"] or any(char.isspace() for char in hit["subject"])
            or not 0 <= hit["identity"] <= 100
            or not 0 <= hit["query_coverage"] <= 100
            or hit["evalue"] < 0 or hit["bitscore"] < 0
            or hit["query_length"] <= 0 or hit["alignment_length"] <= 0
            or not 0 <= hit["mismatches"] <= hit["alignment_length"]
            or not 0 <= hit["gap_opens"] <= hit["alignment_length"]
            or not 1 <= hit["query_start"] <= hit["query_length"]
            or not 1 <= hit["query_end"] <= hit["query_length"]
            or hit["subject_start"] <= 0 or hit["subject_end"] <= 0
            or query_span > hit["alignment_length"]
            or subject_span > hit["alignment_length"]):
        raise BlastError(f"BLAST 第 {line_number} 行含不合法的比对数值或坐标。")
    # Full-query exact matches must have no gaps and span every query base;
    # rounded coverage or a 100%-identical short local match alone is insufficient.
    hit["full_length_exact"] = (
        hit["identity"] == 100.0 and hit["mismatches"] == 0
        and hit["gap_opens"] == 0
        and query_span == subject_span == hit["alignment_length"] == hit["query_length"]
    )
    return hit


def parse_blast_results(path, expected_query_ids):
    """Read a complete outfmt-7 report produced with ``BLAST_OUTFMT``.

    Headers and the final processed-query count distinguish genuine zero hits
    from missing or truncated output. Data are assigned using their ``qseqid``
    field, never the user's ID or a description in a comment. No caller-owned
    results are mutated on a parse failure.
    """
    expected = list(expected_query_ids)
    if (not expected or any(not isinstance(value, str) or not value
                            or any(c.isspace() for c in value) for value in expected)
            or len(set(expected)) != len(expected)):
        raise ValueError("预期的 BLAST 查询 ID 必须是非空、互不重复且不含空白的文本。")
    results = {query_id: [] for query_id in expected}
    seen = set()
    current_query = None
    declared_count = None
    processed_count = None
    expected_lengths = {}

    def finish_query():
        if current_query is not None:
            if declared_count is None:
                raise BlastError(f"查询 {current_query} 缺少命中数，BLAST 输出可能不完整。")
            actual = len(results[current_query])
            if actual != declared_count:
                raise BlastError(
                    f"查询 {current_query} 声明 {declared_count} 条比对，"
                    f"实际读取 {actual} 条，BLAST 输出可能已截断。"
                )

    try:
        with open(path, "r", encoding="utf-8-sig") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                header = re.match(r"^#\s*Query:\s*(\S+)", line, re.IGNORECASE)
                hits_header = re.fullmatch(r"#\s*(\d+)\s+hits?\s+found\s*", line, re.IGNORECASE)
                footer = re.fullmatch(
                    r"#\s*BLAST\s+processed\s+(\d+)\s+quer(?:y|ies)\s*", line, re.IGNORECASE
                )
                if header:
                    if processed_count is not None:
                        raise BlastError("BLAST 完成标记后仍有查询，输出不是一次完整检索。")
                    finish_query()
                    current_query = _query_id(header.group(1), results, line_number)
                    if current_query in seen:
                        raise BlastError(f"BLAST 输出重复出现查询 {current_query}。")
                    seen.add(current_query)
                    declared_count = None
                elif hits_header:
                    if current_query is None or declared_count is not None:
                        raise BlastError(f"BLAST 第 {line_number} 行的命中数没有对应的唯一查询。")
                    declared_count = int(hits_header.group(1))
                elif footer:
                    if processed_count is not None:
                        raise BlastError("BLAST 输出包含重复的完成标记。")
                    finish_query()
                    processed_count = int(footer.group(1))
                    current_query = None
                elif line.startswith("#"):
                    continue  # Version, database, and field-label comments vary by BLAST release.
                else:
                    if current_query is None or declared_count is None:
                        raise BlastError(f"BLAST 第 {line_number} 行的比对缺少查询或命中数标题。")
                    parts = line.split("\t")
                    hit = _parse_hit(parts, line_number)
                    query_id = _query_id(parts[0], results, line_number)
                    if query_id != current_query:
                        raise BlastError(f"BLAST 第 {line_number} 行的查询 ID 与所属查询不一致。")
                    if query_id in expected_lengths and expected_lengths[query_id] != hit["query_length"]:
                        raise BlastError(f"BLAST 查询 {query_id} 的长度在输出中不一致。")
                    expected_lengths[query_id] = hit["query_length"]
                    results[query_id].append(hit)
    except (OSError, UnicodeError) as exc:
        raise BlastError(f"无法读取 BLAST 结果：{exc}") from exc
    finish_query()
    if processed_count is None:
        raise BlastError("BLAST 输出为空或缺少完成标记，不能按“无匹配”处理。")
    if seen != set(expected) or processed_count != len(expected):
        raise BlastError(
            f"BLAST 预期 {len(expected)} 条查询，实际读取 {len(seen)} 条，"
            f"完成标记为 {processed_count} 条；输出不完整。"
        )
    return results


def summarize_hits(hits, max_target_seqs=DEFAULT_MAX_TARGET_SEQS):
    """Describe observed alignments without inferring biological specificity."""
    if (not isinstance(max_target_seqs, Integral)
            or isinstance(max_target_seqs, bool) or max_target_seqs < 1):
        raise ValueError("BLAST 报告的最大目标序列数必须是正整数。")
    hits = list(hits)
    subjects = {hit["subject"] for hit in hits}
    full_matches = {hit["subject"] for hit in hits if hit.get("full_length_exact", False)}
    count = len(subjects)
    if not count:
        label = "无匹配（当前参数）"
    elif count == 1:
        label = "单一序列命中（待核实）"
    else:
        label = "多序列命中（待核实）"
    limit_reached = count >= max_target_seqs
    if limit_reached:
        label += "；达到报告上限"
    return {
        "specificity": label,
        "blast_hits_count": count,
        "blast_hsp_count": len(hits),
        "blast_full_match_count": len(full_matches),
        "blast_hit_limit_reached": limit_reached,
        "blast_report_limit": max_target_seqs,
    }
