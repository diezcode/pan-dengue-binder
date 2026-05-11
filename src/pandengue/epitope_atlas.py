"""Target-atlas generation and epitope hypothesis ranking."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Sequence

from .config import get_nested, resolve_path
from .metadata import sha256_file
from .sequences import require_sequence_manifest


CANDIDATE_COLUMNS = (
    "rank",
    "epitope_id",
    "target_class",
    "e_start",
    "e_end",
    "e_positions",
    "polyprotein_positions",
    "consensus_residues",
    "residue_count",
    "support_count_min",
    "support_count_mean",
    "entropy_mean",
    "gap_fraction_mean",
    "consensus_frequency_mean",
    "structure_evidence",
    "visualizable_residues",
    "conservation_score",
    "structure_score",
    "functional_literature_score",
    "ade_context_score",
    "designability_score",
    "template_score",
    "total_score",
    "advance_eligible",
    "rationale",
    "failure_modes",
)


@dataclass(frozen=True)
class ConservedResidue:
    alignment_column: int
    e_position: int
    polyprotein_position: int | None
    consensus_residue: str
    consensus_frequency: float
    entropy: float
    gap_fraction: float
    support_count: int
    support_accessions: tuple[str, ...]
    serotype_variability: str = ""


@dataclass(frozen=True)
class StructureAnnotation:
    e_position: int
    chain_id: str = ""
    residue_id: str = ""
    template_id: str = ""
    context: str = ""
    exposed: bool | None = None
    relative_sasa: float | None = None
    glycan_distance: float | None = None
    x: float | None = None
    y: float | None = None
    z: float | None = None

    @property
    def has_coordinates(self) -> bool:
        return self.x is not None and self.y is not None and self.z is not None


@dataclass(frozen=True)
class EpitopeCandidate:
    epitope_id: str
    positions: tuple[ConservedResidue, ...]
    target_class: str
    score_components: dict[str, float]
    total_score: float
    advance_eligible: bool
    rationale: tuple[str, ...]
    failure_modes: tuple[str, ...]
    visualizable_residues: tuple[str, ...]
    structure_evidence: str


def plan_target_atlas(config: dict) -> dict:
    outputs_dir = _target_atlas_dir(config)
    return {
        "action": "rank_target_hypotheses_from_conservation_and_structure",
        "dry_run_safe": True,
        "requires": {
            "manifest": str(_sequence_manifest_path_from_config(config)),
            "conservation": str(_conservation_path(config)),
            "optional_structure_mapping": str(_structure_mapping_path(config)),
        },
        "outputs": [
            str(outputs_dir / "epitope_candidates.tsv"),
            str(outputs_dir / "ranked_epitopes.json"),
            str(resolve_path(config, "paths.reports") / "target_atlas.md"),
        ],
        "decision_gate": "advance_2_to_4_epitopes_only_with_sequence_structure_and_literature_support",
    }


def build_target_atlas(config: dict) -> dict:
    """Generate ranked epitope hypotheses from conservation and optional structure evidence."""
    manifest_path = require_sequence_manifest(config)
    conservation_path = _conservation_path(config)
    if not conservation_path.exists():
        raise FileNotFoundError(
            "Phase 4 requires Phase 2 conservation output before target ranking: "
            f"{conservation_path}"
        )

    structure_mapping_path = _structure_mapping_path(config)
    target_atlas_dir = _target_atlas_dir(config)
    reports_dir = resolve_path(config, "paths.reports")
    target_atlas_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    conservation_rows = read_conservation_rows(conservation_path)
    structure_annotations = (
        read_structure_mapping(structure_mapping_path)
        if structure_mapping_path.exists()
        else {}
    )
    candidates = generate_epitope_candidates(conservation_rows, structure_annotations, config)
    ranked_candidates = rank_epitope_candidates(candidates)
    selected = select_pilot_epitopes(ranked_candidates, config)

    candidates_path = target_atlas_dir / str(
        get_nested(config, "target_atlas.candidates_output", "epitope_candidates.tsv")
    )
    ranked_json_path = target_atlas_dir / str(
        get_nested(config, "target_atlas.ranked_output", "ranked_epitopes.json")
    )
    report_path = reports_dir / str(
        get_nested(config, "target_atlas.report_output", "target_atlas.md")
    )

    write_candidates_tsv(candidates_path, ranked_candidates)
    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest_path),
        "conservation": str(conservation_path),
        "structure_mapping": str(structure_mapping_path) if structure_mapping_path.exists() else "",
        "structure_mapping_present": structure_mapping_path.exists(),
        "candidate_count": len(ranked_candidates),
        "selected_count": len(selected),
        "selected_epitope_ids": [candidate.epitope_id for candidate in selected],
        "decision_gate": decision_gate_summary(ranked_candidates, selected, config),
        "candidates": [
            candidate_to_dict(candidate, rank=index + 1)
            for index, candidate in enumerate(ranked_candidates)
        ],
        "outputs": {
            "epitope_candidates": str(candidates_path),
            "ranked_epitopes": str(ranked_json_path),
            "report": str(report_path),
        },
        "input_checksums": _target_atlas_input_checksums(
            manifest_path,
            conservation_path,
            structure_mapping_path,
        ),
    }
    ranked_json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    report_path.write_text(target_atlas_report(payload), encoding="utf-8")
    return payload


def read_conservation_rows(path: str | Path) -> list[ConservedResidue]:
    residues: list[ConservedResidue] = []
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            e_position = _optional_int(row.get("reference_e_position", ""))
            if e_position is None:
                continue
            residues.append(
                ConservedResidue(
                    alignment_column=int(row.get("alignment_column", "0") or 0),
                    e_position=e_position,
                    polyprotein_position=_optional_int(row.get("reference_polyprotein_position", "")),
                    consensus_residue=row.get("consensus_residue", ""),
                    consensus_frequency=_float(row.get("consensus_frequency", "0")),
                    entropy=_float(row.get("shannon_entropy", "0")),
                    gap_fraction=_float(row.get("gap_fraction", "0")),
                    support_count=int(float(row.get("support_count", "0") or 0)),
                    support_accessions=tuple(
                        accession
                        for accession in row.get("support_accessions", "").split(";")
                        if accession
                    ),
                    serotype_variability=row.get("serotype_variability", ""),
                )
            )
    return residues


def read_structure_mapping(path: str | Path) -> dict[int, StructureAnnotation]:
    annotations: dict[int, StructureAnnotation] = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            e_position = _optional_int(
                row.get("reference_e_position", "")
                or row.get("e_position", "")
                or row.get("e_protein_position", "")
            )
            if e_position is None:
                continue
            annotations[e_position] = StructureAnnotation(
                e_position=e_position,
                chain_id=(
                    row.get("chain_id", "")
                    or row.get("auth_asym_id", "")
                    or row.get("label_asym_id", "")
                ),
                residue_id=(
                    row.get("residue_id", "")
                    or row.get("auth_seq_id", "")
                    or row.get("label_seq_id", "")
                ),
                template_id=row.get("template_id", "") or row.get("pdb_id", ""),
                context=row.get("context", "") or row.get("biological_context", ""),
                exposed=_optional_bool(
                    row.get("surface_accessible", "")
                    or row.get("exposed", "")
                    or row.get("is_exposed", "")
                    or row.get("mature_virion_accessible", "")
                    or row.get("recombinant_accessible", "")
                ),
                relative_sasa=_optional_float(
                    row.get("relative_sasa", "")
                    or row.get("sasa_relative", "")
                    or row.get("solvent_accessibility", "")
                ),
                glycan_distance=_optional_float(
                    row.get("glycan_distance_angstrom", "")
                    or row.get("distance_to_glycan", "")
                ),
                x=_optional_float(row.get("x", "")),
                y=_optional_float(row.get("y", "")),
                z=_optional_float(row.get("z", "")),
            )
    return annotations


def generate_epitope_candidates(
    conservation_rows: Sequence[ConservedResidue],
    structure_annotations: dict[int, StructureAnnotation],
    config: dict,
) -> list[EpitopeCandidate]:
    eligible = [
        residue
        for residue in conservation_rows
        if residue.support_count >= int(get_nested(config, "target_atlas.min_support_count", 1))
        and residue.entropy <= float(get_nested(config, "target_atlas.max_entropy", 0.35))
        and residue.gap_fraction <= float(get_nested(config, "target_atlas.max_gap_fraction", 0.10))
        and residue.consensus_frequency >= float(
            get_nested(config, "target_atlas.min_consensus_frequency", 0.85)
        )
    ]
    clusters = cluster_residues(eligible, structure_annotations, config)
    minimum_patch_size = int(get_nested(config, "target_atlas.min_patch_size", 3))

    candidates: list[EpitopeCandidate] = []
    for index, cluster in enumerate(clusters, start=1):
        if len(cluster) < minimum_patch_size:
            continue
        target_class = classify_target_patch(cluster, config)
        candidate_id = f"epitope_{index:03d}"
        candidates.append(
            score_epitope_candidate(
                candidate_id,
                cluster,
                target_class,
                structure_annotations,
                config,
            )
        )
    return candidates


def cluster_residues(
    residues: Sequence[ConservedResidue],
    structure_annotations: dict[int, StructureAnnotation],
    config: dict,
) -> list[tuple[ConservedResidue, ...]]:
    if not residues:
        return []
    ordered = sorted(residues, key=lambda residue: residue.e_position)
    max_sequence_gap = int(get_nested(config, "target_atlas.max_patch_gap", 3))
    distance_threshold = float(
        get_nested(config, "target_atlas.structure_cluster_distance_angstrom", 12.0)
    )

    parent = {residue.e_position: residue.e_position for residue in ordered}

    def find(position: int) -> int:
        while parent[position] != position:
            parent[position] = parent[parent[position]]
            position = parent[position]
        return position

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left, right in zip(ordered, ordered[1:]):
        if right.e_position - left.e_position <= max_sequence_gap:
            union(left.e_position, right.e_position)

    coordinate_residues = [
        residue
        for residue in ordered
        if structure_annotations.get(residue.e_position)
        and structure_annotations[residue.e_position].has_coordinates
    ]
    for index, left in enumerate(coordinate_residues):
        for right in coordinate_residues[index + 1 :]:
            distance = _structure_distance(
                structure_annotations[left.e_position],
                structure_annotations[right.e_position],
            )
            if distance is not None and distance <= distance_threshold:
                union(left.e_position, right.e_position)

    grouped: dict[int, list[ConservedResidue]] = defaultdict(list)
    for residue in ordered:
        grouped[find(residue.e_position)].append(residue)
    return [
        tuple(sorted(group, key=lambda residue: residue.e_position))
        for group in sorted(grouped.values(), key=lambda group: min(residue.e_position for residue in group))
    ]


def classify_target_patch(cluster: Sequence[ConservedResidue], config: dict) -> str:
    start = min(residue.e_position for residue in cluster)
    end = max(residue.e_position for residue in cluster)
    midpoint = (start + end) / 2

    fusion_start = int(get_nested(config, "target_atlas.fusion_loop_start", 98))
    fusion_end = int(get_nested(config, "target_atlas.fusion_loop_end", 112))
    domain_iii_start = int(get_nested(config, "target_atlas.domain_iii_start", 295))
    domain_iii_end = int(get_nested(config, "target_atlas.domain_iii_end", 395))
    stem_start = int(get_nested(config, "target_atlas.stem_start", 396))

    if _overlaps(start, end, fusion_start, fusion_end):
        return "fusion-loop / E Domain II surface"
    if domain_iii_start <= midpoint <= domain_iii_end:
        return "E Domain III hypothesis"
    if midpoint >= stem_start:
        return "E stem/proximal hypothesis"
    if _overlaps(start, end, 250, 295) or _overlaps(start, end, 145, 180):
        return "E dimer-interface hypothesis"
    return "conserved E-protein surface hypothesis"


def score_epitope_candidate(
    epitope_id: str,
    cluster: Sequence[ConservedResidue],
    target_class: str,
    structure_annotations: dict[int, StructureAnnotation],
    config: dict,
) -> EpitopeCandidate:
    entropy_values = [residue.entropy for residue in cluster]
    gap_values = [residue.gap_fraction for residue in cluster]
    consensus_values = [residue.consensus_frequency for residue in cluster]
    support_values = [residue.support_count for residue in cluster]

    max_entropy = max(float(get_nested(config, "target_atlas.max_entropy", 0.35)), 0.001)
    conservation_score = _clamp(
        (
            _mean(consensus_values)
            + (1.0 - _mean(gap_values))
            + (1.0 - min(_mean(entropy_values) / max_entropy, 1.0))
        )
        / 3
    )
    structure_score = _structure_score(cluster, structure_annotations, config)
    functional_literature_score = _class_prior(target_class, "literature")
    ade_context_score = _class_prior(target_class, "ade")
    designability_score = _designability_score(cluster, target_class, structure_score)
    template_score = _template_score(cluster, structure_annotations)

    components = {
        "conservation": conservation_score,
        "structure": structure_score,
        "functional_literature": functional_literature_score,
        "ade_context": ade_context_score,
        "designability": designability_score,
        "template": template_score,
    }
    total_score = _weighted_total(components, config)
    require_structure = bool(get_nested(config, "target_atlas.require_structure_for_advance", True))
    minimum_score = float(get_nested(config, "target_atlas.min_advance_score", 0.65))
    minimum_structure_score = float(
        get_nested(config, "target_atlas.min_structure_score_for_advance", 0.55)
    )
    advance_eligible = (
        total_score >= minimum_score
        and (not require_structure or structure_score >= minimum_structure_score)
    )

    annotations = [
        structure_annotations[residue.e_position]
        for residue in cluster
        if residue.e_position in structure_annotations
    ]
    visualizable_residues = tuple(
        _visualizable_residue(residue, structure_annotations.get(residue.e_position))
        for residue in cluster
    )
    structure_evidence = _structure_evidence_text(cluster, annotations)
    rationale = _candidate_rationale(
        cluster,
        target_class,
        conservation_score,
        structure_score,
        functional_literature_score,
    )
    failure_modes = _candidate_failure_modes(
        target_class,
        structure_annotations_present=bool(structure_annotations),
        structure_score=structure_score,
        annotations=annotations,
    )

    return EpitopeCandidate(
        epitope_id=epitope_id,
        positions=tuple(cluster),
        target_class=target_class,
        score_components=components,
        total_score=total_score,
        advance_eligible=advance_eligible,
        rationale=tuple(rationale),
        failure_modes=tuple(failure_modes),
        visualizable_residues=visualizable_residues,
        structure_evidence=structure_evidence,
    )


def rank_epitope_candidates(candidates: Sequence[EpitopeCandidate]) -> list[EpitopeCandidate]:
    return sorted(
        candidates,
        key=lambda candidate: (
            -candidate.total_score,
            min(residue.e_position for residue in candidate.positions),
            candidate.epitope_id,
        ),
    )


def select_pilot_epitopes(candidates: Sequence[EpitopeCandidate], config: dict) -> list[EpitopeCandidate]:
    advance_min = int(get_nested(config, "target_atlas.advance_min", 2))
    advance_max = int(get_nested(config, "target_atlas.advance_max", 4))
    eligible = [candidate for candidate in candidates if candidate.advance_eligible]
    if len(eligible) < advance_min:
        return []
    return eligible[:advance_max]


def decision_gate_summary(
    candidates: Sequence[EpitopeCandidate],
    selected: Sequence[EpitopeCandidate],
    config: dict,
) -> dict[str, object]:
    advance_min = int(get_nested(config, "target_atlas.advance_min", 2))
    advance_max = int(get_nested(config, "target_atlas.advance_max", 4))
    eligible_count = sum(1 for candidate in candidates if candidate.advance_eligible)
    if selected:
        message = f"advance {len(selected)} epitope hypotheses for pilot design"
        status = "advance"
    elif eligible_count == 1:
        message = (
            "hold: only one epitope has enough support, so pan-serotype robustness should not be claimed"
        )
        status = "hold"
    else:
        message = "hold: expand sequence, structure, or literature evidence before pilot design"
        status = "hold"
    return {
        "status": status,
        "message": message,
        "eligible_count": eligible_count,
        "advance_min": advance_min,
        "advance_max": advance_max,
    }


def write_candidates_tsv(path: str | Path, candidates: Sequence[EpitopeCandidate]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CANDIDATE_COLUMNS), delimiter="\t")
        writer.writeheader()
        for index, candidate in enumerate(candidates, start=1):
            writer.writerow(candidate_to_tsv_row(candidate, rank=index))
    return output_path


def candidate_to_tsv_row(candidate: EpitopeCandidate, *, rank: int) -> dict[str, object]:
    positions = candidate.positions
    e_positions = [residue.e_position for residue in positions]
    polyprotein_positions = [
        str(residue.polyprotein_position)
        for residue in positions
        if residue.polyprotein_position is not None
    ]
    return {
        "rank": rank,
        "epitope_id": candidate.epitope_id,
        "target_class": candidate.target_class,
        "e_start": min(e_positions),
        "e_end": max(e_positions),
        "e_positions": ",".join(str(position) for position in e_positions),
        "polyprotein_positions": ",".join(polyprotein_positions),
        "consensus_residues": "".join(residue.consensus_residue for residue in positions),
        "residue_count": len(positions),
        "support_count_min": min(residue.support_count for residue in positions),
        "support_count_mean": _format_float(_mean([residue.support_count for residue in positions])),
        "entropy_mean": _format_float(_mean([residue.entropy for residue in positions])),
        "gap_fraction_mean": _format_float(_mean([residue.gap_fraction for residue in positions])),
        "consensus_frequency_mean": _format_float(
            _mean([residue.consensus_frequency for residue in positions])
        ),
        "structure_evidence": candidate.structure_evidence,
        "visualizable_residues": ";".join(candidate.visualizable_residues),
        "conservation_score": _format_float(candidate.score_components["conservation"]),
        "structure_score": _format_float(candidate.score_components["structure"]),
        "functional_literature_score": _format_float(
            candidate.score_components["functional_literature"]
        ),
        "ade_context_score": _format_float(candidate.score_components["ade_context"]),
        "designability_score": _format_float(candidate.score_components["designability"]),
        "template_score": _format_float(candidate.score_components["template"]),
        "total_score": _format_float(candidate.total_score),
        "advance_eligible": str(candidate.advance_eligible).lower(),
        "rationale": " ".join(candidate.rationale),
        "failure_modes": " ".join(candidate.failure_modes),
    }


def candidate_to_dict(candidate: EpitopeCandidate, *, rank: int) -> dict[str, object]:
    row = candidate_to_tsv_row(candidate, rank=rank)
    row["score_components"] = candidate.score_components
    row["rationale"] = list(candidate.rationale)
    row["failure_modes"] = list(candidate.failure_modes)
    row["visualizable_residues"] = list(candidate.visualizable_residues)
    return row


def target_atlas_report(payload: dict) -> str:
    decision_gate = payload["decision_gate"]
    candidates = payload["candidates"]
    selected_ids = set(payload["selected_epitope_ids"])
    lines = [
        "# Target Atlas",
        "",
        "Scope: computational epitope hypothesis ranking for non-Fc binder design. This report does not claim efficacy, neutralization, or clinical utility.",
        "",
        "## Decision Gate",
        "",
        f"- Status: `{decision_gate['status']}`",
        f"- Recommendation: {decision_gate['message']}",
        f"- Eligible epitopes: {decision_gate['eligible_count']}",
        f"- Required pilot range: {decision_gate['advance_min']}-{decision_gate['advance_max']}",
        "",
        "## Ranked Epitopes",
        "",
    ]
    if not candidates:
        lines.append("No candidate patches passed the configured conservation filters.")
    else:
        lines.extend(
            [
                "| Rank | Epitope | Class | E positions | Score | Pilot |",
                "|---:|---|---|---|---:|---|",
            ]
        )
        for candidate in candidates:
            selected = "yes" if candidate["epitope_id"] in selected_ids else "no"
            lines.append(
                f"| {candidate['rank']} | `{candidate['epitope_id']}` | {candidate['target_class']} | "
                f"{candidate['e_positions']} | {candidate['total_score']} | {selected} |"
            )

    lines.extend(["", "## Rationale", ""])
    for candidate in candidates[:10]:
        lines.extend(
            [
                f"### {candidate['epitope_id']} - {candidate['target_class']}",
                "",
                f"- E positions: `{candidate['e_positions']}`",
                f"- Visualizable residues: `{';'.join(candidate['visualizable_residues'])}`",
                f"- Score components: conservation `{candidate['conservation_score']}`, structure `{candidate['structure_score']}`, literature/context `{candidate['functional_literature_score']}`, ADE context `{candidate['ade_context_score']}`, designability `{candidate['designability_score']}`, template `{candidate['template_score']}`.",
                f"- Rationale: {' '.join(candidate['rationale'])}",
                f"- Failure modes: {' '.join(candidate['failure_modes'])}",
                "",
            ]
        )

    lines.extend(
        [
            "## Failure Modes Across The Atlas",
            "",
            "- Conservation alone is not evidence that a residue is exposed on a mature virion.",
            "- Recombinant E-protein accessibility may differ from E-dimer or virion-like accessibility.",
            "- Fusion-loop-like targets can be conserved but may carry weak cross-reactive or ADE-adjacent antibody-context risk; non-Fc binder formats remain the default assumption.",
            "- Domain III is treated as one target hypothesis among several, not as a hardcoded endpoint.",
            "",
            "## Outputs",
            "",
        ]
    )
    for label, path in payload["outputs"].items():
        lines.append(f"- `{label}`: `{path}`")
    lines.append("")
    return "\n".join(lines)


def _structure_score(
    cluster: Sequence[ConservedResidue],
    annotations_by_position: dict[int, StructureAnnotation],
    config: dict,
) -> float:
    if not annotations_by_position:
        return float(get_nested(config, "target_atlas.missing_structure_score", 0.20))
    annotations = [
        annotations_by_position.get(residue.e_position)
        for residue in cluster
        if residue.e_position in annotations_by_position
    ]
    if not annotations:
        return float(get_nested(config, "target_atlas.unmapped_structure_score", 0.15))

    exposure_scores: list[float] = []
    for annotation in annotations:
        if annotation.exposed is not None:
            exposure_scores.append(1.0 if annotation.exposed else 0.0)
        elif annotation.relative_sasa is not None:
            exposure_scores.append(_clamp(annotation.relative_sasa))
        else:
            exposure_scores.append(0.5)
    coverage = len(annotations) / len(cluster)
    return _clamp((0.65 * _mean(exposure_scores)) + (0.35 * coverage))


def _template_score(
    cluster: Sequence[ConservedResidue],
    annotations_by_position: dict[int, StructureAnnotation],
) -> float:
    annotations = [
        annotations_by_position.get(residue.e_position)
        for residue in cluster
        if residue.e_position in annotations_by_position
    ]
    if not annotations:
        return 0.20
    mapped_fraction = len(annotations) / len(cluster)
    template_fraction = sum(1 for annotation in annotations if annotation.template_id) / len(cluster)
    chain_fraction = sum(1 for annotation in annotations if annotation.chain_id and annotation.residue_id) / len(cluster)
    return _clamp((mapped_fraction + template_fraction + chain_fraction) / 3)


def _designability_score(
    cluster: Sequence[ConservedResidue],
    target_class: str,
    structure_score: float,
) -> float:
    size = len(cluster)
    if 6 <= size <= 25:
        size_score = 1.0
    elif 3 <= size <= 40:
        size_score = 0.75
    else:
        size_score = 0.45
    class_score = _class_prior(target_class, "designability")
    return _clamp((0.45 * size_score) + (0.35 * structure_score) + (0.20 * class_score))


def _weighted_total(components: dict[str, float], config: dict) -> float:
    weights = {
        "conservation": float(get_nested(config, "target_atlas.weights.conservation", 0.35)),
        "structure": float(get_nested(config, "target_atlas.weights.structure", 0.20)),
        "functional_literature": float(
            get_nested(config, "target_atlas.weights.functional_literature", 0.15)
        ),
        "ade_context": float(get_nested(config, "target_atlas.weights.ade_context", 0.10)),
        "designability": float(get_nested(config, "target_atlas.weights.designability", 0.10)),
        "template": float(get_nested(config, "target_atlas.weights.template", 0.10)),
    }
    total_weight = sum(weights.values()) or 1.0
    return _clamp(sum(components[key] * weights[key] for key in weights) / total_weight)


def _class_prior(target_class: str, category: str) -> float:
    priors = {
        "E Domain III hypothesis": {
            "literature": 0.75,
            "ade": 0.70,
            "designability": 0.80,
        },
        "fusion-loop / E Domain II surface": {
            "literature": 0.55,
            "ade": 0.40,
            "designability": 0.55,
        },
        "E dimer-interface hypothesis": {
            "literature": 0.65,
            "ade": 0.70,
            "designability": 0.65,
        },
        "E stem/proximal hypothesis": {
            "literature": 0.45,
            "ade": 0.60,
            "designability": 0.45,
        },
        "conserved E-protein surface hypothesis": {
            "literature": 0.45,
            "ade": 0.60,
            "designability": 0.60,
        },
    }
    return priors.get(target_class, priors["conserved E-protein surface hypothesis"])[category]


def _candidate_rationale(
    cluster: Sequence[ConservedResidue],
    target_class: str,
    conservation_score: float,
    structure_score: float,
    literature_score: float,
) -> list[str]:
    start = min(residue.e_position for residue in cluster)
    end = max(residue.e_position for residue in cluster)
    return [
        f"Patch spans E positions {start}-{end} and is classified as {target_class}.",
        f"Mean conservation support is strong enough for ranking (conservation score {conservation_score:.2f}).",
        f"Structure evidence score is {structure_score:.2f}; literature/context prior is {literature_score:.2f}.",
    ]


def _candidate_failure_modes(
    target_class: str,
    *,
    structure_annotations_present: bool,
    structure_score: float,
    annotations: Sequence[StructureAnnotation],
) -> list[str]:
    failures: list[str] = []
    if not structure_annotations_present:
        failures.append("No structure mapping was available, so mature-virion exposure is unresolved.")
    elif structure_score < 0.55:
        failures.append("Structure evidence is weak or incomplete for this patch.")
    if "fusion-loop" in target_class:
        failures.append("Fusion-loop conservation can coincide with cross-reactive weak-binding/ADE-context risk.")
    if "stem" in target_class:
        failures.append("Stem/proximal accessibility may depend strongly on virion state and breathing.")
    if any(
        annotation.glycan_distance is not None and annotation.glycan_distance < 6.0
        for annotation in annotations
    ):
        failures.append("Nearby glycan annotation may reduce binder access or complicate modeling.")
    if not failures:
        failures.append("Computational ranking still requires later structural and experimental validation.")
    return failures


def _structure_evidence_text(
    cluster: Sequence[ConservedResidue],
    annotations: Sequence[StructureAnnotation],
) -> str:
    if not annotations:
        return "no_structure_mapping"
    exposed = sum(1 for annotation in annotations if annotation.exposed is True)
    mapped = len(annotations)
    templates = sorted({annotation.template_id for annotation in annotations if annotation.template_id})
    contexts = sorted({annotation.context for annotation in annotations if annotation.context})
    parts = [f"mapped={mapped}/{len(cluster)}", f"exposed={exposed}/{mapped}"]
    if templates:
        parts.append(f"templates={','.join(templates)}")
    if contexts:
        parts.append(f"contexts={','.join(contexts)}")
    return ";".join(parts)


def _visualizable_residue(
    residue: ConservedResidue,
    annotation: StructureAnnotation | None,
) -> str:
    if annotation and annotation.chain_id and annotation.residue_id:
        template = f"{annotation.template_id}:" if annotation.template_id else ""
        return f"{template}{annotation.chain_id}:{annotation.residue_id}"
    return f"E:{residue.e_position}"


def _target_atlas_input_checksums(
    manifest_path: Path,
    conservation_path: Path,
    structure_mapping_path: Path,
) -> dict[str, str]:
    checksums = {
        "manifest": sha256_file(manifest_path),
        "conservation": sha256_file(conservation_path),
    }
    if structure_mapping_path.exists():
        checksums["structure_mapping"] = sha256_file(structure_mapping_path)
    return checksums


def _structure_distance(
    left: StructureAnnotation,
    right: StructureAnnotation,
) -> float | None:
    if not left.has_coordinates or not right.has_coordinates:
        return None
    assert left.x is not None and left.y is not None and left.z is not None
    assert right.x is not None and right.y is not None and right.z is not None
    return math.sqrt(
        (left.x - right.x) ** 2
        + (left.y - right.y) ** 2
        + (left.z - right.z) ** 2
    )


def _sequence_manifest_path_from_config(config: dict) -> Path:
    processed_sequences = resolve_path(config, "paths.processed_sequences")
    return processed_sequences / str(
        get_nested(config, "sequence_data.manifest_output", "sequence_manifest.tsv")
    )


def _conservation_path(config: dict) -> Path:
    return _target_atlas_dir(config) / str(
        get_nested(config, "alignment.conservation_output", "conservation.tsv")
    )


def _structure_mapping_path(config: dict) -> Path:
    return _target_atlas_dir(config) / str(
        get_nested(config, "target_atlas.structure_mapping_input", "structure_mapping.tsv")
    )


def _target_atlas_dir(config: dict) -> Path:
    return resolve_path(config, "paths.outputs") / "target_atlas"


def _overlaps(left_start: int, left_end: int, right_start: int, right_end: int) -> bool:
    return left_start <= right_end and right_start <= left_end


def _mean(values: Sequence[float | int]) -> float:
    return sum(float(value) for value in values) / len(values) if values else 0.0


def _clamp(value: float, lower: float = 0.0, upper: float = 1.0) -> float:
    return min(max(value, lower), upper)


def _format_float(value: float) -> str:
    return f"{value:.6f}"


def _float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def _optional_float(value: str) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _optional_int(value: str) -> int | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def _optional_bool(value: str) -> bool | None:
    if value is None or str(value).strip() == "":
        return None
    normalized = str(value).strip().lower()
    if normalized in {"true", "yes", "y", "1", "exposed", "accessible"}:
        return True
    if normalized in {"false", "no", "n", "0", "buried", "inaccessible"}:
        return False
    return None
