from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import sys


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pandengue.cli import run_stage
from pandengue.sequences import (
    SequenceGateError,
    SourceSequence,
    build_sequence_dataset,
    extract_e_protein,
    read_fasta,
    read_manifest_tsv,
)


def phase1_config(tmp_path: Path, fasta_path: Path, manifest_path: Path | None = None) -> dict:
    return {
        "project": {"name": "phase1-test"},
        "paths": {
            "sequences": str(tmp_path / "raw" / "sequences"),
            "processed_sequences": str(tmp_path / "processed" / "sequences"),
            "reports": str(tmp_path / "reports"),
            "outputs": str(tmp_path / "outputs"),
        },
        "data_sources": {"uniprot_reference_polyproteins": {}},
        "pipeline": {
            "dry_run_default": True,
            "minimum_sequences_per_serotype": 1,
            "run_metadata_filename": "run_metadata.json",
        },
        "sequence_data": {
            "include_uniprot_references": False,
            "local_fasta": str(fasta_path),
            "local_manifest": str(manifest_path or ""),
            "cache_subdir": "uniprot",
            "polyprotein_output": "denv_polyproteins.fasta",
            "e_protein_output": "denv_e_proteins.fasta",
            "manifest_output": "sequence_manifest.tsv",
            "summary_output": "sequence_summary.json",
            "e_protein": {
                "start_1_based": 2,
                "end_1_based_inclusive": 7,
            },
            "validation": {
                "min_e_protein_length": 4,
                "max_e_protein_length": 12,
            },
        },
    }


class Phase1SequenceLayerTests(unittest.TestCase):
    def test_build_sequence_dataset_deduplicates_and_reports_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "bulk.fasta"
            fasta_path.write_text(
                ">ACC1 DENV-1 envelope isolate\nACDEFG\n"
                ">ACC1_DUP DENV-1 envelope duplicate\nACDEFG\n"
                ">ACC2 DENV-2 envelope isolate\nCDEFGH\n"
                ">ACC3 DENV-3 envelope isolate\nDEFGHI\n"
                ">ACC4 DENV-4 envelope isolate\nEFGHIK\n",
                encoding="utf-8",
            )
            config = phase1_config(tmp_path, fasta_path)

            summary = build_sequence_dataset(config, fetch_remote=False)

            self.assertTrue(summary["gates"]["passed"])
            self.assertEqual(summary["counts"]["included_records_total"], 5)
            self.assertEqual(summary["counts"]["unique_e_sequences_total"], 4)
            self.assertEqual(summary["counts"]["unique_e_sequences_by_serotype"]["DENV-1"], 1)

            e_records = read_fasta(summary["outputs"]["e_protein_fasta"])
            self.assertEqual(len(e_records), 4)
            self.assertIn("duplicate_count=2", e_records[0].description)

            manifest_rows = read_manifest_tsv(summary["outputs"]["manifest"])
            self.assertEqual(len(manifest_rows), 5)
            self.assertTrue(all(row["accession"] for row in manifest_rows))
            self.assertTrue(all(row["source"] for row in manifest_rows))

            report_text = Path(summary["outputs"]["report"]).read_text(encoding="utf-8")
            self.assertIn("| DENV-4 | 1 | 1 | 1 |", report_text)

    def test_annotation_coordinates_are_preferred_for_e_extraction(self) -> None:
        source = SourceSequence(
            accession="ANN1",
            source="local_manifest",
            serotype="DENV-1",
            description="DENV-1 polyprotein",
            sequence="MABCDEQRST",
            protein="polyprotein",
            metadata={"e_start_1_based": "2", "e_end_1_based": "6"},
        )
        extracted = extract_e_protein(source, phase1_config(Path("tmp"), Path("unused.fasta")))
        self.assertEqual(extracted.sequence, "ABCDE")
        self.assertEqual(extracted.method, "annotation_coordinates")
        self.assertEqual(extracted.start_1_based, 2)
        self.assertEqual(extracted.end_1_based, 6)

    def test_sequence_gates_require_all_four_serotypes(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "incomplete.fasta"
            fasta_path.write_text(">ACC1 DENV-1 envelope isolate\nACDEFG\n", encoding="utf-8")
            config = phase1_config(tmp_path, fasta_path)

            with self.assertRaises(SequenceGateError) as context:
                build_sequence_dataset(config, fetch_remote=False)

            self.assertFalse(context.exception.summary["gates"]["passed"])
            self.assertIn("DENV-4 has 0", context.exception.summary["gates"]["message"])
            self.assertTrue(Path(context.exception.summary["outputs"]["report"]).exists())

    def test_cli_fetch_sequences_executes_with_local_input(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "bulk.fasta"
            fasta_path.write_text(
                ">A1 DENV-1\nAAAA\n"
                ">A2 DENV-2\nCCCC\n"
                ">A3 DENV-3\nDDDD\n"
                ">A4 DENV-4\nEEEE\n",
                encoding="utf-8",
            )
            result = run_stage("fetch-sequences", phase1_config(tmp_path, fasta_path), dry_run=False)

            self.assertEqual(result["status"], "completed")
            self.assertTrue(Path(result["execution"]["outputs"]["manifest"]).exists())
            self.assertTrue(Path(result["metadata_path"]).exists())

    def test_alignment_stage_refuses_to_run_without_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "bulk.fasta"
            fasta_path.write_text(">A1 DENV-1\nAAAA\n", encoding="utf-8")
            config = phase1_config(tmp_path, fasta_path)

            with self.assertRaises(FileNotFoundError):
                run_stage("build-alignment", config, dry_run=False)


if __name__ == "__main__":
    unittest.main()
