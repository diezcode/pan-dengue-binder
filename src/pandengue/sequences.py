"""Sequence data ingestion, normalization, and Phase 1 dataset gates."""

from __future__ import annotations

import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Sequence
from urllib.error import URLError
from urllib.request import urlopen

from .config import get_nested, repo_root, resolve_path
from .metadata import sha256_file


SEROTYPES = ("DENV-1", "DENV-2", "DENV-3", "DENV-4")

MANIFEST_COLUMNS = (
    "normalized_id",
    "accession",
    "source",
    "serotype",
    "protein",
    "sequence_role",
    "description",
    "host",
    "country_region",
    "collection_date",
    "retrieval_date",
    "raw_sequence_length",
    "normalized_sequence_length",
    "extraction_method",
    "region_start_1_based",
    "region_end_1_based",
    "deduplicated_sequence_id",
    "duplicate_count",
    "filter_status",
    "filter_reason",
)


class SequenceGateError(ValueError):
    """Raised when a sequence dataset does not satisfy configured gates."""

    def __init__(self, message: str, summary: dict | None = None) -> None:
        super().__init__(message)
        self.summary = summary or {}


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class SourceSequence:
    accession: str
    source: str
    serotype: str
    description: str
    sequence: str
    protein: str = ""
    sequence_role: str = "bulk_isolate"
    host: str = ""
    country_region: str = ""
    collection_date: str = ""
    retrieval_date: str = ""
    metadata: dict[str, str] | None = None


@dataclass(frozen=True)
class ExtractedSequence:
    source: SourceSequence
    sequence: str
    method: str
    start_1_based: int
    end_1_based: int


def parse_fasta(text: str) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    sequence_lines: list[str] = []

    def flush() -> None:
        if header is None:
            return
        sequence = "".join(sequence_lines).replace(" ", "").upper()
        identifier = header.split()[0]
        records.append(FastaRecord(identifier, header, sequence))

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].strip()
            sequence_lines = []
        else:
            if header is None:
                raise ValueError("FASTA sequence data appeared before a header")
            sequence_lines.append(line)
    flush()
    return records


def read_fasta(path: str | Path) -> list[FastaRecord]:
    return parse_fasta(Path(path).read_text(encoding="utf-8"))


def write_fasta(path: str | Path, records: Iterable[FastaRecord], *, line_width: int = 80) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for record in records:
        lines.append(f">{record.description}")
        sequence = record.sequence
        lines.extend(sequence[index : index + line_width] for index in range(0, len(sequence), line_width))
    output_path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    return output_path


def parse_manifest_tsv(text: str) -> list[dict[str, str]]:
    lines = [line.rstrip("\n") for line in text.splitlines() if line.strip()]
    if not lines:
        return []
    headers = lines[0].split("\t")
    rows: list[dict[str, str]] = []
    for line in lines[1:]:
        values = line.split("\t")
        row = {header: values[index] if index < len(values) else "" for index, header in enumerate(headers)}
        rows.append(row)
    return rows


def read_manifest_tsv(path: str | Path) -> list[dict[str, str]]:
    return parse_manifest_tsv(Path(path).read_text(encoding="utf-8"))


def write_manifest_tsv(path: str | Path, rows: Iterable[dict[str, object]]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rendered = ["\t".join(MANIFEST_COLUMNS)]
    for row in rows:
        rendered.append(
            "\t".join(_clean_tsv_value(str(row.get(column, ""))) for column in MANIFEST_COLUMNS)
        )
    output_path.write_text("\n".join(rendered).rstrip() + "\n", encoding="utf-8")
    return output_path


def conserved_pattern(records: Iterable[FastaRecord]) -> str:
    sequences = [record.sequence for record in records]
    if not sequences:
        return ""
    lengths = {len(sequence) for sequence in sequences}
    if len(lengths) != 1:
        raise ValueError("Conserved pattern requires equal-length aligned sequences")

    pattern: list[str] = []
    for residues in zip(*sequences):
        first = residues[0]
        pattern.append(first if all(residue == first for residue in residues) else "-")
    return "".join(pattern)


def find_shared_windows(
    reference: str,
    comparison_sequences: Iterable[str],
    *,
    min_length: int = 10,
) -> list[str]:
    if min_length < 1:
        raise ValueError("min_length must be positive")

    comparisons = list(comparison_sequences)
    windows: list[str] = []
    index = 0
    while index <= len(reference) - min_length:
        window = reference[index : index + min_length]
        if all(window in sequence for sequence in comparisons):
            end = index + min_length
            while end < len(reference):
                candidate = reference[index : end + 1]
                if all(candidate in sequence for sequence in comparisons):
                    end += 1
                    continue
                break
            windows.append(reference[index:end])
            index = end
        else:
            index += 1
    return windows


def plan_sequence_fetch(config: dict) -> dict:
    references = config.get("data_sources", {}).get("uniprot_reference_polyproteins", {})
    processed_dir = str(get_nested(config, "paths.processed_sequences", "data/processed/sequences"))
    return {
        "action": "fetch_cache_and_normalize_sequence_dataset",
        "dry_run_safe": True,
        "reference_count": len(references),
        "references": references,
        "supports": [
            "curated UniProt reference polyproteins",
            "local bulk isolate FASTA files",
            "local sequence metadata manifests",
            "E-protein extraction and identical-sequence deduplication",
        ],
        "outputs": [
            "data/raw/sequences/",
            f"{processed_dir}/denv_e_proteins.fasta",
            f"{processed_dir}/sequence_manifest.tsv",
            f"{processed_dir}/sequence_summary.json",
            "reports/sequence_data.md",
        ],
    }


def sequence_manifest_path(config: dict) -> Path:
    return _processed_sequences_dir(config) / str(
        get_nested(config, "sequence_data.manifest_output", "sequence_manifest.tsv")
    )


def require_sequence_manifest(config: dict) -> Path:
    manifest_path = sequence_manifest_path(config)
    if not manifest_path.exists():
        raise FileNotFoundError(
            "Sequence manifest is required before conservation or target-atlas outputs can be produced: "
            f"{manifest_path}"
        )
    rows = read_manifest_tsv(manifest_path)
    usable_rows = [
        row
        for row in rows
        if row.get("filter_status") == "included"
        and row.get("deduplicated_sequence_id")
        and row.get("accession")
        and row.get("source")
    ]
    if not usable_rows:
        raise ValueError(f"Sequence manifest contains no included E-protein records: {manifest_path}")
    return manifest_path


def build_sequence_dataset(config: dict, *, fetch_remote: bool = True) -> dict:
    """Build the Phase 1 normalized E-protein FASTA, manifest, and report."""
    raw_sequences_dir = _raw_sequences_dir(config)
    processed_dir = _processed_sequences_dir(config)
    reports_dir = resolve_path(config, "paths.reports")
    raw_sequences_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    retrieval_date = datetime.now(timezone.utc).date().isoformat()
    manifest_rows = _read_configured_manifest_rows(config)
    source_records = list(_load_source_sequences(config, manifest_rows, retrieval_date, fetch_remote))

    extracted: list[ExtractedSequence] = []
    output_rows: list[dict[str, object]] = []
    filter_counter: Counter[str] = Counter()

    for source_record in source_records:
        try:
            e_protein = extract_e_protein(source_record, config)
            _validate_e_protein(e_protein.sequence, config)
        except ValueError as exc:
            filter_counter[str(exc)] += 1
            output_rows.append(_manifest_row_for_excluded(source_record, str(exc)))
            continue
        extracted.append(e_protein)

    deduped_records, dedup_index = _deduplicate_extracted_sequences(extracted)
    for e_protein in extracted:
        dedup_id, duplicate_count = dedup_index[(e_protein.source.serotype, e_protein.sequence)]
        output_rows.append(_manifest_row_for_included(e_protein, dedup_id, duplicate_count))

    polyprotein_path = processed_dir / str(
        get_nested(config, "sequence_data.polyprotein_output", "denv_polyproteins.fasta")
    )
    e_protein_path = processed_dir / str(
        get_nested(config, "sequence_data.e_protein_output", "denv_e_proteins.fasta")
    )
    manifest_path = sequence_manifest_path(config)
    summary_path = processed_dir / str(
        get_nested(config, "sequence_data.summary_output", "sequence_summary.json")
    )
    report_path = reports_dir / "sequence_data.md"

    write_fasta(
        polyprotein_path,
        _source_records_as_fasta(source_records),
    )
    write_fasta(e_protein_path, deduped_records)
    write_manifest_tsv(manifest_path, output_rows)

    summary = _build_summary(
        config=config,
        source_records=source_records,
        extracted=extracted,
        deduped_records=deduped_records,
        filter_counter=filter_counter,
        outputs={
            "polyprotein_fasta": str(polyprotein_path),
            "e_protein_fasta": str(e_protein_path),
            "manifest": str(manifest_path),
            "summary": str(summary_path),
            "report": str(report_path),
        },
    )
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(sequence_data_report(summary), encoding="utf-8")

    if not summary["gates"]["passed"]:
        raise SequenceGateError(summary["gates"]["message"], summary)
    return summary


def extract_e_protein(source: SourceSequence, config: dict) -> ExtractedSequence:
    """Extract a mature E-protein sequence from a source record."""
    sequence = _clean_sequence(source.sequence)
    metadata = source.metadata or {}

    annotation_coords = _coordinate_pair(
        metadata,
        (
            ("e_start_1_based", "e_end_1_based"),
            ("e_start", "e_end"),
            ("protein_start_1_based", "protein_end_1_based"),
            ("region_start_1_based", "region_end_1_based"),
            ("start_1_based", "end_1_based"),
        ),
    )
    if annotation_coords is not None:
        start, end = annotation_coords
        return _slice_sequence(
            source,
            sequence,
            start,
            end,
            "annotation_coordinates",
        )

    if _looks_like_e_protein(source.protein, source.description) and _within_e_length(sequence, config):
        return ExtractedSequence(
            source=source,
            sequence=sequence,
            method="provided_e_protein",
            start_1_based=1,
            end_1_based=len(sequence),
        )

    configured_start = int(get_nested(config, "sequence_data.e_protein.start_1_based", 281))
    configured_end = int(get_nested(config, "sequence_data.e_protein.end_1_based_inclusive", 775))
    if len(sequence) >= configured_end:
        return _slice_sequence(
            source,
            sequence,
            configured_start,
            configured_end,
            "reference_guided_config_coordinates",
        )

    if _within_e_length(sequence, config):
        return ExtractedSequence(
            source=source,
            sequence=sequence,
            method="length_based_e_protein",
            start_1_based=1,
            end_1_based=len(sequence),
        )

    raise ValueError("missing_e_protein_coordinates")


def infer_serotype(*values: str) -> str:
    text = " ".join(value for value in values if value)
    patterns = (
        r"\bDENV[\s_-]*([1-4])\b",
        r"\bdengue(?:\s+virus)?(?:\s+type)?[\s_-]*([1-4])\b",
        r"\bserotype[\s_:=/-]*([1-4])\b",
        r"\btype[\s_:=/-]*([1-4])\b",
    )
    for pattern in patterns:
        match = re.search(pattern, text, flags=re.IGNORECASE)
        if match:
            return f"DENV-{match.group(1)}"
    return ""


def accession_from_identifier(identifier: str) -> str:
    parts = [part for part in identifier.split("|") if part]
    if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
        return parts[1]
    return parts[0] if parts else identifier


def sequence_data_report(summary: dict) -> str:
    lines = [
        "# Sequence Data Report",
        "",
        "Scope: computational sequence provenance and filtering for target-atlas analysis.",
        "",
        "## Dataset Counts",
        "",
        "| Serotype | Source records | Included records | Unique E sequences |",
        "|---|---:|---:|---:|",
    ]
    for serotype in SEROTYPES:
        source_count = summary["counts"]["source_records_by_serotype"].get(serotype, 0)
        included_count = summary["counts"]["included_records_by_serotype"].get(serotype, 0)
        unique_count = summary["counts"]["unique_e_sequences_by_serotype"].get(serotype, 0)
        lines.append(f"| {serotype} | {source_count} | {included_count} | {unique_count} |")

    lines.extend(["", "## Filters", ""])
    if summary["filters"]:
        for reason, count in sorted(summary["filters"].items()):
            lines.append(f"- `{reason}`: {count}")
    else:
        lines.append("- No source records were excluded.")

    lines.extend(
        [
            "",
            "## Dataset Gates",
            "",
            f"- Status: {'passed' if summary['gates']['passed'] else 'failed'}",
            f"- Minimum unique E sequences per serotype: {summary['gates']['minimum_unique_e_sequences_per_serotype']}",
            f"- Message: {summary['gates']['message']}",
            "",
            "## Outputs",
            "",
        ]
    )
    for label, path in summary["outputs"].items():
        lines.append(f"- `{label}`: `{path}`")

    lines.extend(
        [
            "",
            "This report does not make therapeutic, clinical, or experimental claims.",
            "",
        ]
    )
    return "\n".join(lines)


def _load_source_sequences(
    config: dict,
    manifest_rows: Sequence[dict[str, str]],
    retrieval_date: str,
    fetch_remote: bool,
) -> Iterator[SourceSequence]:
    metadata_by_key = _metadata_lookup(manifest_rows)
    local_paths = [
        ("sequence_data.local_fasta", "local_fasta", "bulk_isolate"),
        ("sequence_data.local_bulk_fasta", "local_fasta", "bulk_isolate"),
        ("sequence_data.local_curated_fasta", "local_curated_fasta", "curated_reference"),
    ]
    for dotted_key, default_source, sequence_role in local_paths:
        local_path = _optional_path(config, dotted_key)
        if local_path is None:
            continue
        for record in read_fasta(local_path):
            yield _source_sequence_from_fasta(
                record,
                metadata_by_key,
                default_source=default_source,
                default_role=sequence_role,
                default_retrieval_date=retrieval_date,
            )

    if bool(get_nested(config, "sequence_data.include_uniprot_references", True)):
        references = config.get("data_sources", {}).get("uniprot_reference_polyproteins", {})
        cache_dir = _raw_sequences_dir(config) / str(
            get_nested(config, "sequence_data.cache_subdir", "uniprot")
        )
        cache_dir.mkdir(parents=True, exist_ok=True)
        for serotype, accession in references.items():
            fasta_text = _read_or_fetch_uniprot_fasta(
                accession=str(accession),
                cache_dir=cache_dir,
                fetch_remote=fetch_remote,
            )
            for record in parse_fasta(fasta_text):
                metadata = {
                    "accession": str(accession),
                    "source": "UniProt",
                    "serotype": str(serotype),
                    "protein": "polyprotein",
                    "retrieval_date": retrieval_date,
                }
                yield _source_sequence_from_fasta(
                    record,
                    metadata_by_key,
                    default_source="UniProt",
                    default_role="curated_reference",
                    default_retrieval_date=retrieval_date,
                    overrides=metadata,
                )


def _read_or_fetch_uniprot_fasta(accession: str, cache_dir: Path, fetch_remote: bool) -> str:
    cache_path = cache_dir / f"{accession}.fasta"
    if cache_path.exists():
        return cache_path.read_text(encoding="utf-8")
    if not fetch_remote:
        raise FileNotFoundError(f"Cached UniProt FASTA not found: {cache_path}")

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    try:
        with urlopen(url, timeout=30) as response:
            text = response.read().decode("utf-8")
    except URLError as exc:
        raise RuntimeError(f"Unable to fetch UniProt FASTA for {accession}: {exc}") from exc
    if not text.startswith(">"):
        raise RuntimeError(f"UniProt FASTA response for {accession} did not contain FASTA data")
    cache_path.write_text(text, encoding="utf-8")
    return text


def _source_sequence_from_fasta(
    record: FastaRecord,
    metadata_by_key: dict[str, dict[str, str]],
    *,
    default_source: str,
    default_role: str,
    default_retrieval_date: str,
    overrides: dict[str, str] | None = None,
) -> SourceSequence:
    accession = accession_from_identifier(record.identifier)
    row = {}
    for key in (accession, record.identifier):
        if key in metadata_by_key:
            row = metadata_by_key[key]
            break
    if overrides:
        row = {**row, **overrides}

    description = row.get("description") or record.description
    serotype = _normalize_serotype(
        row.get("serotype", "") or infer_serotype(record.identifier, record.description, description)
    )
    source = row.get("source") or row.get("source_database") or default_source
    return SourceSequence(
        accession=row.get("accession") or accession,
        source=source,
        serotype=serotype,
        description=description,
        sequence=_clean_sequence(record.sequence),
        protein=row.get("protein") or row.get("product", ""),
        sequence_role=row.get("sequence_role") or default_role,
        host=row.get("host", ""),
        country_region=row.get("country_region") or row.get("country") or row.get("region", ""),
        collection_date=row.get("collection_date", ""),
        retrieval_date=row.get("retrieval_date") or default_retrieval_date,
        metadata=row,
    )


def _deduplicate_extracted_sequences(
    extracted: Sequence[ExtractedSequence],
) -> tuple[list[FastaRecord], dict[tuple[str, str], tuple[str, int]]]:
    grouped: dict[tuple[str, str], list[ExtractedSequence]] = defaultdict(list)
    for e_protein in extracted:
        grouped[(e_protein.source.serotype, e_protein.sequence)].append(e_protein)

    by_serotype_counter: Counter[str] = Counter()
    fasta_records: list[FastaRecord] = []
    dedup_index: dict[tuple[str, str], tuple[str, int]] = {}
    for serotype, sequence in sorted(grouped):
        group = grouped[(serotype, sequence)]
        by_serotype_counter[serotype] += 1
        dedup_id = f"{serotype}|E|{by_serotype_counter[serotype]:04d}"
        duplicate_count = len(group)
        first = group[0].source
        description = (
            f"{dedup_id} accession={first.accession} serotype={serotype} "
            f"source={first.source} duplicate_count={duplicate_count}"
        )
        fasta_records.append(FastaRecord(dedup_id, description, sequence))
        dedup_index[(serotype, sequence)] = (dedup_id, duplicate_count)
    return fasta_records, dedup_index


def _build_summary(
    *,
    config: dict,
    source_records: Sequence[SourceSequence],
    extracted: Sequence[ExtractedSequence],
    deduped_records: Sequence[FastaRecord],
    filter_counter: Counter[str],
    outputs: dict[str, str],
) -> dict:
    source_counts = Counter(record.serotype or "unknown" for record in source_records)
    included_counts = Counter(record.source.serotype for record in extracted)
    unique_counts = Counter(record.identifier.split("|", 1)[0] for record in deduped_records)
    minimum = int(get_nested(config, "pipeline.minimum_sequences_per_serotype", 1))
    gate_messages: list[str] = []
    for serotype in SEROTYPES:
        if unique_counts.get(serotype, 0) < minimum:
            gate_messages.append(
                f"{serotype} has {unique_counts.get(serotype, 0)} unique E sequences; requires {minimum}"
            )
    passed = not gate_messages
    return {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "counts": {
            "source_records_total": len(source_records),
            "included_records_total": len(extracted),
            "unique_e_sequences_total": len(deduped_records),
            "source_records_by_serotype": dict(sorted(source_counts.items())),
            "included_records_by_serotype": dict(sorted(included_counts.items())),
            "unique_e_sequences_by_serotype": dict(sorted(unique_counts.items())),
        },
        "filters": dict(sorted(filter_counter.items())),
        "gates": {
            "passed": passed,
            "minimum_unique_e_sequences_per_serotype": minimum,
            "required_serotypes": list(SEROTYPES),
            "message": "passed" if passed else "; ".join(gate_messages),
        },
        "outputs": outputs,
        "input_checksums": _existing_output_checksums(outputs),
    }


def _existing_output_checksums(outputs: dict[str, str]) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for label, raw_path in outputs.items():
        path = Path(raw_path)
        if path.exists() and path.is_file() and label != "summary":
            checksums[label] = sha256_file(path)
    return checksums


def _manifest_row_for_included(
    e_protein: ExtractedSequence,
    dedup_id: str,
    duplicate_count: int,
) -> dict[str, object]:
    source = e_protein.source
    return {
        "normalized_id": f"{source.accession}|E",
        "accession": source.accession,
        "source": source.source,
        "serotype": source.serotype,
        "protein": "E",
        "sequence_role": source.sequence_role,
        "description": source.description,
        "host": source.host,
        "country_region": source.country_region,
        "collection_date": source.collection_date,
        "retrieval_date": source.retrieval_date,
        "raw_sequence_length": len(source.sequence),
        "normalized_sequence_length": len(e_protein.sequence),
        "extraction_method": e_protein.method,
        "region_start_1_based": e_protein.start_1_based,
        "region_end_1_based": e_protein.end_1_based,
        "deduplicated_sequence_id": dedup_id,
        "duplicate_count": duplicate_count,
        "filter_status": "included",
        "filter_reason": "",
    }


def _manifest_row_for_excluded(source: SourceSequence, reason: str) -> dict[str, object]:
    return {
        "normalized_id": "",
        "accession": source.accession,
        "source": source.source,
        "serotype": source.serotype,
        "protein": source.protein,
        "sequence_role": source.sequence_role,
        "description": source.description,
        "host": source.host,
        "country_region": source.country_region,
        "collection_date": source.collection_date,
        "retrieval_date": source.retrieval_date,
        "raw_sequence_length": len(source.sequence),
        "normalized_sequence_length": 0,
        "extraction_method": "",
        "region_start_1_based": "",
        "region_end_1_based": "",
        "deduplicated_sequence_id": "",
        "duplicate_count": 0,
        "filter_status": "excluded",
        "filter_reason": reason,
    }


def _source_records_as_fasta(source_records: Iterable[SourceSequence]) -> Iterator[FastaRecord]:
    for record in source_records:
        description = (
            f"{record.accession} serotype={record.serotype or 'unknown'} "
            f"source={record.source} role={record.sequence_role} {record.description}"
        )
        yield FastaRecord(record.accession, description, record.sequence)


def _validate_e_protein(sequence: str, config: dict) -> None:
    if not _within_e_length(sequence, config):
        raise ValueError("e_protein_length_out_of_range")
    invalid = sorted(set(sequence) - set("ABCDEFGHIKLMNPQRSTVWXYZ*-UO"))
    if invalid:
        raise ValueError("invalid_amino_acid_characters")


def _within_e_length(sequence: str, config: dict) -> bool:
    minimum = int(get_nested(config, "sequence_data.validation.min_e_protein_length", 400))
    maximum = int(get_nested(config, "sequence_data.validation.max_e_protein_length", 600))
    return minimum <= len(sequence) <= maximum


def _slice_sequence(
    source: SourceSequence,
    sequence: str,
    start_1_based: int,
    end_1_based: int,
    method: str,
) -> ExtractedSequence:
    if start_1_based < 1 or end_1_based < start_1_based or end_1_based > len(sequence):
        raise ValueError("invalid_e_protein_coordinates")
    return ExtractedSequence(
        source=source,
        sequence=sequence[start_1_based - 1 : end_1_based],
        method=method,
        start_1_based=start_1_based,
        end_1_based=end_1_based,
    )


def _coordinate_pair(rows: dict[str, str], key_pairs: Sequence[tuple[str, str]]) -> tuple[int, int] | None:
    for start_key, end_key in key_pairs:
        raw_start = rows.get(start_key, "")
        raw_end = rows.get(end_key, "")
        if raw_start and raw_end:
            try:
                return int(raw_start), int(raw_end)
            except ValueError as exc:
                raise ValueError("invalid_e_protein_coordinates") from exc
    return None


def _looks_like_e_protein(protein: str, description: str) -> bool:
    text = f"{protein} {description}".lower()
    return any(
        marker in text
        for marker in (
            "envelope protein",
            "envelope glycoprotein",
            "protein e",
            " e protein",
            " e glycoprotein",
        )
    )


def _metadata_lookup(rows: Sequence[dict[str, str]]) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for row in rows:
        normalized = {key.strip(): value.strip() for key, value in row.items()}
        for key in (
            normalized.get("accession", ""),
            normalized.get("identifier", ""),
            normalized.get("normalized_id", ""),
        ):
            if key:
                lookup[key] = normalized
    return lookup


def _read_configured_manifest_rows(config: dict) -> list[dict[str, str]]:
    manifest_path = _optional_path(config, "sequence_data.local_manifest")
    if manifest_path is None:
        return []
    return read_manifest_tsv(manifest_path)


def _normalize_serotype(value: str) -> str:
    if value in SEROTYPES:
        return value
    inferred = infer_serotype(value)
    return inferred if inferred in SEROTYPES else ""


def _clean_sequence(sequence: str) -> str:
    return re.sub(r"\s+", "", sequence).upper()


def _clean_tsv_value(value: str) -> str:
    return value.replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def _optional_path(config: dict, dotted_key: str) -> Path | None:
    raw = get_nested(config, dotted_key, "")
    if raw is None or str(raw).strip() == "":
        return None
    path = Path(str(raw))
    if path.is_absolute():
        return path
    return repo_root() / path


def _raw_sequences_dir(config: dict) -> Path:
    return resolve_path(config, "paths.sequences")


def _processed_sequences_dir(config: dict) -> Path:
    return resolve_path(config, "paths.processed_sequences")
