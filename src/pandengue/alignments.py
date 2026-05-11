"""Multiple sequence alignment and conservation-table generation."""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

from .config import get_nested, repo_root, resolve_path
from .metadata import sha256_file
from .sequences import (
    FastaRecord,
    read_fasta,
    read_manifest_tsv,
    require_sequence_manifest,
    sequence_manifest_path,
    write_fasta,
)


CONSERVATION_COLUMNS = (
    "alignment_column",
    "reference_sequence_id",
    "reference_serotype",
    "reference_residue",
    "reference_e_position",
    "reference_polyprotein_position",
    "consensus_residue",
    "consensus_frequency",
    "shannon_entropy",
    "gap_fraction",
    "support_count",
    "global_residue_frequencies",
    "serotype_variability",
    "support_accessions",
)

SEROTYPE_CONSERVATION_COLUMNS = (
    "alignment_column",
    "serotype",
    "reference_e_position",
    "reference_polyprotein_position",
    "consensus_residue",
    "consensus_frequency",
    "shannon_entropy",
    "gap_fraction",
    "support_count",
    "unique_residue_count",
    "residue_frequencies",
    "support_accessions",
)


@dataclass(frozen=True)
class AlignmentRecordMetadata:
    record_id: str
    serotype: str
    accessions: tuple[str, ...]
    source_count: int
    region_start_1_based: int | None = None
    region_end_1_based: int | None = None


def plan_alignment(config: dict) -> dict:
    tools = config.get("external_tools", {})
    manifest_path = sequence_manifest_path(config)
    return {
        "action": "build_multiple_sequence_alignment_and_conservation_tables",
        "dry_run_safe": True,
        "requires_manifest": str(manifest_path),
        "manifest_exists": manifest_path.exists(),
        "preferred_tools": {
            "mafft": tools.get("mafft", {}).get("executable", "mafft"),
            "clustal_omega": tools.get("clustal_omega", {}).get("executable", "clustalo"),
        },
        "fallback": "reference_guided_python_msa",
        "outputs": [
            str(_processed_alignments_dir(config)),
            str(_target_atlas_dir(config) / "conservation.tsv"),
            str(_target_atlas_dir(config) / "conservation_by_serotype.tsv"),
            str(_target_atlas_dir(config) / "alignment_summary.json"),
        ],
    }


def build_alignment_dataset(config: dict) -> dict:
    """Build an MSA and export Phase 2 conservation metrics."""
    manifest_path = require_sequence_manifest(config)
    input_fasta = _e_protein_fasta_path(config)
    if not input_fasta.exists():
        raise FileNotFoundError(f"Normalized E-protein FASTA is required: {input_fasta}")

    processed_alignments_dir = _processed_alignments_dir(config)
    target_atlas_dir = _target_atlas_dir(config)
    processed_alignments_dir.mkdir(parents=True, exist_ok=True)
    target_atlas_dir.mkdir(parents=True, exist_ok=True)

    alignment_path = processed_alignments_dir / str(
        get_nested(config, "alignment.output_fasta", "denv_e_proteins.aligned.fasta")
    )
    alignment_method = build_alignment(config, input_fasta, alignment_path)
    aligned_records, parser_name = read_alignment(alignment_path)
    manifest_rows = read_manifest_tsv(manifest_path)
    metadata_by_id = alignment_metadata_from_manifest(manifest_rows)

    reference_id = select_reference_record(aligned_records, metadata_by_id, config)
    conservation_rows, serotype_rows, metrics_summary = compute_conservation_tables(
        aligned_records,
        metadata_by_id,
        reference_id=reference_id,
    )

    conservation_path = target_atlas_dir / str(
        get_nested(config, "alignment.conservation_output", "conservation.tsv")
    )
    conservation_by_serotype_path = target_atlas_dir / str(
        get_nested(config, "alignment.conservation_by_serotype_output", "conservation_by_serotype.tsv")
    )
    summary_path = target_atlas_dir / str(
        get_nested(config, "alignment.summary_output", "alignment_summary.json")
    )

    _write_tsv(conservation_path, CONSERVATION_COLUMNS, conservation_rows)
    _write_tsv(
        conservation_by_serotype_path,
        SEROTYPE_CONSERVATION_COLUMNS,
        serotype_rows,
    )

    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "alignment_method": alignment_method,
        "alignment_parser": parser_name,
        "input_fasta": str(input_fasta),
        "manifest": str(manifest_path),
        "alignment_fasta": str(alignment_path),
        "record_count": len(aligned_records),
        "alignment_length": len(aligned_records[0].sequence) if aligned_records else 0,
        "reference_sequence_id": reference_id,
        **metrics_summary,
        "outputs": {
            "alignment_fasta": str(alignment_path),
            "conservation": str(conservation_path),
            "conservation_by_serotype": str(conservation_by_serotype_path),
            "summary": str(summary_path),
        },
        "input_checksums": {
            "input_fasta": sha256_file(input_fasta),
            "manifest": sha256_file(manifest_path),
            "alignment_fasta": sha256_file(alignment_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_alignment(config: dict, input_fasta: str | Path, output_fasta: str | Path) -> str:
    """Run the configured aligner, falling back to the pure-Python aligner when allowed."""
    method = str(get_nested(config, "alignment.method", "auto")).lower()
    input_path = Path(input_fasta)
    output_path = Path(output_fasta)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if method in {"mafft", "auto"}:
        executable = str(get_nested(config, "external_tools.mafft.executable", "mafft"))
        if method == "mafft" or shutil.which(executable):
            if shutil.which(executable):
                _run_mafft(executable, input_path, output_path)
                return "mafft"
            raise FileNotFoundError(f"MAFFT executable not found: {executable}")

    if method in {"clustal_omega", "clustalo", "auto"}:
        executable = str(get_nested(config, "external_tools.clustal_omega.executable", "clustalo"))
        if method in {"clustal_omega", "clustalo"} or shutil.which(executable):
            if shutil.which(executable):
                _run_clustal_omega(executable, input_path, output_path)
                return "clustal_omega"
            raise FileNotFoundError(f"Clustal Omega executable not found: {executable}")

    if method in {"python", "fallback", "auto"}:
        records = read_fasta(input_path)
        limit = int(get_nested(config, "alignment.python_fallback_max_records", 100))
        if len(records) > limit:
            raise RuntimeError(
                f"Python fallback MSA is limited to {limit} records; configure MAFFT or Clustal Omega"
            )
        reference_id = str(get_nested(config, "alignment.reference_sequence_id", ""))
        aligned_records = reference_guided_msa(records, reference_id=reference_id or None)
        write_fasta(output_path, aligned_records)
        return "python_reference_guided"

    raise ValueError(f"Unknown alignment method: {method}")


def build_mafft_command(executable: str, input_fasta: str | Path) -> list[str]:
    return [executable, "--auto", str(input_fasta)]


def build_clustal_omega_command(
    executable: str,
    input_fasta: str | Path,
    output_fasta: str | Path,
) -> list[str]:
    return [
        executable,
        "-i",
        str(input_fasta),
        "-o",
        str(output_fasta),
        "--force",
        "--outfmt=fasta",
    ]


def read_alignment(path: str | Path) -> tuple[list[FastaRecord], str]:
    """Parse an aligned FASTA, preferring Biopython AlignIO when installed."""
    alignment_path = Path(path)
    try:
        from Bio import AlignIO  # type: ignore
    except ImportError:
        records = read_fasta(alignment_path)
        _validate_alignment_records(records)
        return records, "standard_library_fasta"

    alignment = AlignIO.read(str(alignment_path), "fasta")
    records = [
        FastaRecord(record.id, record.description, str(record.seq).upper())
        for record in alignment
    ]
    _validate_alignment_records(records)
    return records, "biopython_alignio"


def reference_guided_msa(
    records: Sequence[FastaRecord],
    *,
    reference_id: str | None = None,
) -> list[FastaRecord]:
    """Small deterministic MSA fallback anchored to one reference sequence."""
    if not records:
        return []
    reference = _select_reference_for_python_msa(records, reference_id)
    reference_sequence = reference.sequence.replace("-", "")

    max_insertions = [0] * (len(reference_sequence) + 1)
    per_record_slots: dict[str, tuple[list[str], list[list[str]]]] = {}

    for record in records:
        query_sequence = record.sequence.replace("-", "")
        if record.identifier == reference.identifier:
            aligned_reference, aligned_query = reference_sequence, reference_sequence
        else:
            aligned_reference, aligned_query = needleman_wunsch(reference_sequence, query_sequence)
        residues, insertions = _pairwise_to_reference_slots(
            aligned_reference,
            aligned_query,
            reference_length=len(reference_sequence),
        )
        per_record_slots[record.identifier] = (residues, insertions)
        for slot_index, inserted in enumerate(insertions):
            max_insertions[slot_index] = max(max_insertions[slot_index], len(inserted))

    aligned_records: list[FastaRecord] = []
    for record in records:
        residues, insertions = per_record_slots[record.identifier]
        aligned_sequence_parts: list[str] = []
        aligned_sequence_parts.extend(_padded_insertions(insertions[0], max_insertions[0]))
        for position, residue in enumerate(residues, start=1):
            aligned_sequence_parts.append(residue)
            aligned_sequence_parts.extend(
                _padded_insertions(insertions[position], max_insertions[position])
            )
        aligned_records.append(
            FastaRecord(
                identifier=record.identifier,
                description=record.description,
                sequence="".join(aligned_sequence_parts),
            )
        )

    _validate_alignment_records(aligned_records)
    return aligned_records


def needleman_wunsch(reference: str, query: str) -> tuple[str, str]:
    """Global pairwise alignment used by the fallback MSA."""
    match_score = 2
    mismatch_score = -1
    gap_score = -1
    rows = len(reference) + 1
    cols = len(query) + 1
    scores = [[0] * cols for _ in range(rows)]
    traceback = [[""] * cols for _ in range(rows)]

    for row in range(1, rows):
        scores[row][0] = row * gap_score
        traceback[row][0] = "up"
    for col in range(1, cols):
        scores[0][col] = col * gap_score
        traceback[0][col] = "left"

    for row in range(1, rows):
        for col in range(1, cols):
            diag = scores[row - 1][col - 1] + (
                match_score if reference[row - 1] == query[col - 1] else mismatch_score
            )
            up = scores[row - 1][col] + gap_score
            left = scores[row][col - 1] + gap_score
            best = max(diag, up, left)
            scores[row][col] = best
            if best == diag:
                traceback[row][col] = "diag"
            elif best == up:
                traceback[row][col] = "up"
            else:
                traceback[row][col] = "left"

    aligned_reference: list[str] = []
    aligned_query: list[str] = []
    row = len(reference)
    col = len(query)
    while row > 0 or col > 0:
        direction = traceback[row][col]
        if direction == "diag":
            aligned_reference.append(reference[row - 1])
            aligned_query.append(query[col - 1])
            row -= 1
            col -= 1
        elif direction == "up":
            aligned_reference.append(reference[row - 1])
            aligned_query.append("-")
            row -= 1
        else:
            aligned_reference.append("-")
            aligned_query.append(query[col - 1])
            col -= 1

    return "".join(reversed(aligned_reference)), "".join(reversed(aligned_query))


def alignment_metadata_from_manifest(rows: Iterable[dict[str, str]]) -> dict[str, AlignmentRecordMetadata]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("filter_status") != "included":
            continue
        record_id = row.get("deduplicated_sequence_id", "")
        if record_id:
            grouped[record_id].append(row)

    metadata: dict[str, AlignmentRecordMetadata] = {}
    for record_id, group in grouped.items():
        first = group[0]
        accessions = tuple(row.get("accession", "") for row in group if row.get("accession"))
        metadata[record_id] = AlignmentRecordMetadata(
            record_id=record_id,
            serotype=first.get("serotype", ""),
            accessions=accessions,
            source_count=max(1, len(accessions)),
            region_start_1_based=_optional_int(first.get("region_start_1_based", "")),
            region_end_1_based=_optional_int(first.get("region_end_1_based", "")),
        )
    return metadata


def select_reference_record(
    aligned_records: Sequence[FastaRecord],
    metadata_by_id: dict[str, AlignmentRecordMetadata],
    config: dict,
) -> str:
    configured = str(get_nested(config, "alignment.reference_sequence_id", ""))
    identifiers = {record.identifier for record in aligned_records}
    if configured:
        if configured not in identifiers:
            raise ValueError(f"Configured reference sequence is not in alignment: {configured}")
        return configured

    reference_serotype = str(get_nested(config, "alignment.reference_serotype", "DENV-1"))
    for record in aligned_records:
        metadata = metadata_by_id.get(record.identifier)
        if metadata and metadata.serotype == reference_serotype:
            return record.identifier
    return aligned_records[0].identifier


def compute_conservation_tables(
    aligned_records: Sequence[FastaRecord],
    metadata_by_id: dict[str, AlignmentRecordMetadata],
    *,
    reference_id: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    _validate_alignment_records(aligned_records)
    record_by_id = {record.identifier: record for record in aligned_records}
    if reference_id not in record_by_id:
        raise ValueError(f"Reference record is not present in alignment: {reference_id}")

    reference_record = record_by_id[reference_id]
    reference_metadata = metadata_by_id.get(reference_id)
    reference_serotype = reference_metadata.serotype if reference_metadata else ""
    reference_region_start = (
        reference_metadata.region_start_1_based if reference_metadata else None
    )

    alignment_length = len(reference_record.sequence)
    global_rows: list[dict[str, object]] = []
    serotype_rows: list[dict[str, object]] = []
    reference_e_position = 0
    insertion_columns = 0

    for column_index in range(alignment_length):
        reference_residue = reference_record.sequence[column_index]
        if _is_gap(reference_residue):
            display_e_position: int | str = ""
            display_polyprotein_position: int | str = ""
            insertion_columns += 1
        else:
            reference_e_position += 1
            display_e_position = reference_e_position
            display_polyprotein_position = (
                reference_region_start + reference_e_position - 1
                if reference_region_start is not None
                else ""
            )

        column_entries = [
            (record, record.sequence[column_index], _metadata_for_record(record, metadata_by_id))
            for record in aligned_records
        ]
        global_stats = _weighted_column_stats(column_entries)
        per_serotype_stats: dict[str, dict[str, object]] = {}

        for serotype in sorted({entry[2].serotype for entry in column_entries if entry[2].serotype}):
            serotype_entries = [entry for entry in column_entries if entry[2].serotype == serotype]
            per_serotype_stats[serotype] = _weighted_column_stats(serotype_entries)
            stats = per_serotype_stats[serotype]
            serotype_rows.append(
                {
                    "alignment_column": column_index + 1,
                    "serotype": serotype,
                    "reference_e_position": display_e_position,
                    "reference_polyprotein_position": display_polyprotein_position,
                    "consensus_residue": stats["consensus_residue"],
                    "consensus_frequency": _format_float(stats["consensus_frequency"]),
                    "shannon_entropy": _format_float(stats["shannon_entropy"]),
                    "gap_fraction": _format_float(stats["gap_fraction"]),
                    "support_count": stats["support_count"],
                    "unique_residue_count": stats["unique_residue_count"],
                    "residue_frequencies": _format_frequencies(stats["residue_frequencies"]),
                    "support_accessions": ";".join(stats["support_accessions"]),
                }
            )

        global_rows.append(
            {
                "alignment_column": column_index + 1,
                "reference_sequence_id": reference_id,
                "reference_serotype": reference_serotype,
                "reference_residue": reference_residue,
                "reference_e_position": display_e_position,
                "reference_polyprotein_position": display_polyprotein_position,
                "consensus_residue": global_stats["consensus_residue"],
                "consensus_frequency": _format_float(global_stats["consensus_frequency"]),
                "shannon_entropy": _format_float(global_stats["shannon_entropy"]),
                "gap_fraction": _format_float(global_stats["gap_fraction"]),
                "support_count": global_stats["support_count"],
                "global_residue_frequencies": _format_frequencies(global_stats["residue_frequencies"]),
                "serotype_variability": _format_serotype_variability(per_serotype_stats),
                "support_accessions": ";".join(global_stats["support_accessions"]),
            }
        )

    metadata_counts = Counter(
        _metadata_for_record(record, metadata_by_id).serotype or "unknown"
        for record in aligned_records
    )
    support_counts: Counter[str] = Counter()
    for record in aligned_records:
        metadata = _metadata_for_record(record, metadata_by_id)
        support_counts[metadata.serotype or "unknown"] += metadata.source_count

    summary = {
        "reference_serotype": reference_serotype,
        "reference_region_start_1_based": reference_region_start,
        "reference_columns_numbered": reference_e_position,
        "insertions_relative_to_reference": insertion_columns,
        "serotype_record_counts": dict(sorted(metadata_counts.items())),
        "serotype_support_counts": dict(sorted(support_counts.items())),
        "conservation_columns": len(global_rows),
    }
    return global_rows, serotype_rows, summary


def _weighted_column_stats(
    entries: Sequence[tuple[FastaRecord, str, AlignmentRecordMetadata]],
) -> dict[str, object]:
    total_weight = sum(metadata.source_count for _, _, metadata in entries)
    residue_counts: Counter[str] = Counter()
    gap_count = 0
    support_accessions: list[str] = []

    for _, residue, metadata in entries:
        if _is_gap(residue):
            gap_count += metadata.source_count
            continue
        residue_counts[residue] += metadata.source_count
        support_accessions.extend(metadata.accessions or (metadata.record_id,))

    support_count = sum(residue_counts.values())
    consensus_residue = ""
    consensus_frequency = 0.0
    if residue_counts:
        consensus_residue, consensus_count = sorted(
            residue_counts.items(),
            key=lambda item: (-item[1], item[0]),
        )[0]
        consensus_frequency = consensus_count / support_count

    frequencies = {
        residue: count / support_count
        for residue, count in sorted(residue_counts.items())
    } if support_count else {}

    return {
        "consensus_residue": consensus_residue,
        "consensus_frequency": consensus_frequency,
        "shannon_entropy": _weighted_entropy(residue_counts),
        "gap_fraction": (gap_count / total_weight) if total_weight else 0.0,
        "support_count": support_count,
        "unique_residue_count": len(residue_counts),
        "residue_frequencies": frequencies,
        "support_accessions": sorted(set(support_accessions)),
    }


def _weighted_entropy(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        probability = count / total
        entropy -= probability * math.log2(probability)
    return entropy


def _metadata_for_record(
    record: FastaRecord,
    metadata_by_id: dict[str, AlignmentRecordMetadata],
) -> AlignmentRecordMetadata:
    if record.identifier in metadata_by_id:
        return metadata_by_id[record.identifier]
    serotype = record.identifier.split("|", 1)[0] if record.identifier.startswith("DENV-") else ""
    return AlignmentRecordMetadata(
        record_id=record.identifier,
        serotype=serotype,
        accessions=(record.identifier,),
        source_count=1,
    )


def _run_mafft(executable: str, input_fasta: Path, output_fasta: Path) -> None:
    result = subprocess.run(
        build_mafft_command(executable, input_fasta),
        check=True,
        capture_output=True,
        text=True,
    )
    output_fasta.write_text(result.stdout, encoding="utf-8")


def _run_clustal_omega(executable: str, input_fasta: Path, output_fasta: Path) -> None:
    subprocess.run(
        build_clustal_omega_command(executable, input_fasta, output_fasta),
        check=True,
        capture_output=True,
        text=True,
    )


def _pairwise_to_reference_slots(
    aligned_reference: str,
    aligned_query: str,
    *,
    reference_length: int,
) -> tuple[list[str], list[list[str]]]:
    residues = ["-"] * reference_length
    insertions: list[list[str]] = [[] for _ in range(reference_length + 1)]
    reference_position = 0
    for reference_residue, query_residue in zip(aligned_reference, aligned_query):
        if reference_residue == "-":
            if query_residue != "-":
                insertions[reference_position].append(query_residue)
            continue
        reference_position += 1
        residues[reference_position - 1] = query_residue
    return residues, insertions


def _padded_insertions(insertions: Sequence[str], width: int) -> list[str]:
    return list(insertions) + ["-"] * (width - len(insertions))


def _select_reference_for_python_msa(
    records: Sequence[FastaRecord],
    reference_id: str | None,
) -> FastaRecord:
    if reference_id:
        for record in records:
            if record.identifier == reference_id:
                return record
        raise ValueError(f"Reference sequence is not in input FASTA: {reference_id}")
    for record in records:
        if record.identifier.startswith("DENV-1"):
            return record
    return records[0]


def _validate_alignment_records(records: Sequence[FastaRecord]) -> None:
    if not records:
        raise ValueError("Alignment contains no records")
    lengths = {len(record.sequence) for record in records}
    if len(lengths) != 1:
        raise ValueError("Alignment records must all have the same length")


def _write_tsv(path: str | Path, columns: Sequence[str], rows: Iterable[dict[str, object]]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return output_path


def _format_float(value: object) -> str:
    return f"{float(value):.6f}"


def _format_frequencies(frequencies: object) -> str:
    if not isinstance(frequencies, dict):
        return ""
    return ";".join(f"{residue}:{_format_float(value)}" for residue, value in sorted(frequencies.items()))


def _format_serotype_variability(per_serotype_stats: dict[str, dict[str, object]]) -> str:
    parts: list[str] = []
    for serotype, stats in sorted(per_serotype_stats.items()):
        parts.append(
            f"{serotype}:unique={stats['unique_residue_count']},entropy={_format_float(stats['shannon_entropy'])}"
        )
    return ";".join(parts)


def _optional_int(value: str) -> int | None:
    if not value:
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _is_gap(residue: str) -> bool:
    return residue in {"-", "."}


def _e_protein_fasta_path(config: dict) -> Path:
    processed_sequences = _path_from_config(
        config,
        "paths.processed_sequences",
        Path("data") / "processed" / "sequences",
    )
    return processed_sequences / str(get_nested(config, "sequence_data.e_protein_output", "denv_e_proteins.fasta"))


def _processed_alignments_dir(config: dict) -> Path:
    return _path_from_config(
        config,
        "paths.processed_alignments",
        Path("data") / "processed" / "alignments",
    )


def _target_atlas_dir(config: dict) -> Path:
    outputs_root = resolve_path(config, "paths.outputs")
    return outputs_root / "target_atlas"


def _path_from_config(config: dict, dotted_key: str, default_relative: Path) -> Path:
    raw = get_nested(config, dotted_key, "")
    if raw:
        path = Path(str(raw))
        return path if path.is_absolute() else repo_root() / path
    return repo_root() / default_relative
