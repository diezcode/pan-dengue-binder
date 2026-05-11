from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

import sys


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pandengue.cli import run_stage
from pandengue.epitope_atlas import (
    ConservedResidue,
    StructureAnnotation,
    build_target_atlas,
    cluster_residues,
    generate_epitope_candidates,
)
from pandengue.sequences import MANIFEST_COLUMNS, write_manifest_tsv


CONSERVATION_COLUMNS = [
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
]


def phase4_config(tmp_path: Path) -> dict:
    return {
        "project": {"name": "phase4-test"},
        "paths": {
            "processed_sequences": str(tmp_path / "processed" / "sequences"),
            "outputs": str(tmp_path / "outputs"),
            "reports": str(tmp_path / "reports"),
        },
        "pipeline": {
            "run_metadata_filename": "run_metadata.json",
            "minimum_sequences_per_serotype": 1,
            "dry_run_default": True,
        },
        "sequence_data": {"manifest_output": "sequence_manifest.tsv"},
        "alignment": {"conservation_output": "conservation.tsv"},
        "target_atlas": {
            "structure_mapping_input": "structure_mapping.tsv",
            "candidates_output": "epitope_candidates.tsv",
            "ranked_output": "ranked_epitopes.json",
            "report_output": "target_atlas.md",
            "min_support_count": 2,
            "max_entropy": 0.25,
            "max_gap_fraction": 0.10,
            "min_consensus_frequency": 0.90,
            "min_patch_size": 3,
            "max_patch_gap": 3,
            "structure_cluster_distance_angstrom": 8.0,
            "missing_structure_score": 0.20,
            "unmapped_structure_score": 0.15,
            "min_advance_score": 0.65,
            "min_structure_score_for_advance": 0.55,
            "require_structure_for_advance": True,
            "advance_min": 2,
            "advance_max": 4,
            "fusion_loop_start": 98,
            "fusion_loop_end": 112,
            "domain_iii_start": 295,
            "domain_iii_end": 395,
            "stem_start": 396,
            "weights": {
                "conservation": 0.35,
                "structure": 0.20,
                "functional_literature": 0.15,
                "ade_context": 0.10,
                "designability": 0.10,
                "template": 0.10,
            },
        },
    }


def write_manifest(path: Path) -> None:
    rows = []
    for serotype in ("DENV-1", "DENV-2", "DENV-3", "DENV-4"):
        rows.append(
            {
                "normalized_id": f"{serotype}|ACC|E",
                "accession": f"{serotype}-ACC",
                "source": "unit-test",
                "serotype": serotype,
                "protein": "E",
                "sequence_role": "bulk_isolate",
                "description": f"{serotype} test",
                "host": "",
                "country_region": "",
                "collection_date": "",
                "retrieval_date": "2026-05-11",
                "raw_sequence_length": 5,
                "normalized_sequence_length": 5,
                "extraction_method": "provided_e_protein",
                "region_start_1_based": 281,
                "region_end_1_based": 285,
                "deduplicated_sequence_id": f"{serotype}|E|0001",
                "duplicate_count": 1,
                "filter_status": "included",
                "filter_reason": "",
            }
        )
    write_manifest_tsv(path, rows)


def write_conservation(path: Path, positions: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CONSERVATION_COLUMNS, delimiter="\t")
        writer.writeheader()
        for column, position in enumerate(positions, start=1):
            writer.writerow(
                {
                    "alignment_column": column,
                    "reference_sequence_id": "DENV-1|E|0001",
                    "reference_serotype": "DENV-1",
                    "reference_residue": "A",
                    "reference_e_position": position,
                    "reference_polyprotein_position": 280 + position,
                    "consensus_residue": "A",
                    "consensus_frequency": "0.98",
                    "shannon_entropy": "0.02",
                    "gap_fraction": "0.00",
                    "support_count": "8",
                    "global_residue_frequencies": "A:0.980000",
                    "serotype_variability": "DENV-1:unique=1,entropy=0.000000",
                    "support_accessions": "A1;A2;A3;A4",
                }
            )


def write_structure_mapping(path: Path, positions: list[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "reference_e_position",
            "template_id",
            "chain_id",
            "residue_id",
            "surface_accessible",
            "relative_sasa",
            "context",
            "x",
            "y",
            "z",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for position in positions:
            writer.writerow(
                {
                    "reference_e_position": position,
                    "template_id": "1OAN",
                    "chain_id": "A",
                    "residue_id": str(position),
                    "surface_accessible": "true",
                    "relative_sasa": "0.85",
                    "context": "recombinant_E",
                    "x": str(float(position)),
                    "y": "0.0",
                    "z": "0.0",
                }
            )


class Phase4EpitopeAtlasTests(unittest.TestCase):
    def test_clusters_can_use_structure_distance_for_noncontiguous_residues(self) -> None:
        residues = tuple(
            ConservedResidue(
                alignment_column=index,
                e_position=position,
                polyprotein_position=280 + position,
                consensus_residue="A",
                consensus_frequency=0.98,
                entropy=0.01,
                gap_fraction=0.0,
                support_count=8,
                support_accessions=("A1",),
            )
            for index, position in enumerate((300, 305, 330), start=1)
        )
        annotations = {
            300: StructureAnnotation(300, x=0.0, y=0.0, z=0.0, exposed=True),
            305: StructureAnnotation(305, x=3.0, y=0.0, z=0.0, exposed=True),
            330: StructureAnnotation(330, x=5.0, y=0.0, z=0.0, exposed=True),
        }

        clusters = cluster_residues(residues, annotations, phase4_config(Path("tmp")))

        self.assertEqual(len(clusters), 1)
        self.assertEqual([residue.e_position for residue in clusters[0]], [300, 305, 330])

    def test_candidates_without_structure_are_ranked_but_not_advanced(self) -> None:
        residues = [
            ConservedResidue(
                alignment_column=index,
                e_position=position,
                polyprotein_position=280 + position,
                consensus_residue="A",
                consensus_frequency=0.98,
                entropy=0.02,
                gap_fraction=0.0,
                support_count=8,
                support_accessions=("A1", "A2"),
            )
            for index, position in enumerate((300, 301, 302), start=1)
        ]

        candidates = generate_epitope_candidates(residues, {}, phase4_config(Path("tmp")))

        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0].target_class, "E Domain III hypothesis")
        self.assertFalse(candidates[0].advance_eligible)
        self.assertIn("No structure mapping was available", candidates[0].failure_modes[0])

    def test_build_target_atlas_writes_ranked_outputs_and_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            config = phase4_config(tmp_path)
            manifest_path = tmp_path / "processed" / "sequences" / "sequence_manifest.tsv"
            conservation_path = tmp_path / "outputs" / "target_atlas" / "conservation.tsv"
            structure_path = tmp_path / "outputs" / "target_atlas" / "structure_mapping.tsv"

            write_manifest(manifest_path)
            write_conservation(conservation_path, [145, 146, 147, 300, 301, 302])
            write_structure_mapping(structure_path, [145, 146, 147, 300, 301, 302])

            payload = build_target_atlas(config)

            self.assertEqual(payload["decision_gate"]["status"], "advance")
            self.assertEqual(payload["selected_count"], 2)
            self.assertTrue(Path(payload["outputs"]["epitope_candidates"]).exists())
            self.assertTrue(Path(payload["outputs"]["ranked_epitopes"]).exists())
            report_text = Path(payload["outputs"]["report"]).read_text(encoding="utf-8")
            self.assertIn("Domain III is treated as one target hypothesis", report_text)
            self.assertIn("E dimer-interface hypothesis", report_text)

    def test_cli_build_target_atlas_executes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            config = phase4_config(tmp_path)
            manifest_path = tmp_path / "processed" / "sequences" / "sequence_manifest.tsv"
            conservation_path = tmp_path / "outputs" / "target_atlas" / "conservation.tsv"
            structure_path = tmp_path / "outputs" / "target_atlas" / "structure_mapping.tsv"

            write_manifest(manifest_path)
            write_conservation(conservation_path, [145, 146, 147, 300, 301, 302])
            write_structure_mapping(structure_path, [145, 146, 147, 300, 301, 302])

            result = run_stage("build-target-atlas", config, dry_run=False)

            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["execution"]["selected_count"], 2)
            self.assertTrue(Path(result["metadata_path"]).exists())


if __name__ == "__main__":
    unittest.main()
