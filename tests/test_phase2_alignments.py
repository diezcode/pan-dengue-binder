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

from pandengue.alignments import (
    AlignmentRecordMetadata,
    build_clustal_omega_command,
    build_mafft_command,
    compute_conservation_tables,
    needleman_wunsch,
    reference_guided_msa,
)
from pandengue.cli import run_stage
from pandengue.sequences import FastaRecord, build_sequence_dataset


def phase2_config(tmp_path: Path, fasta_path: Path) -> dict:
    return {
        "project": {"name": "phase2-test"},
        "paths": {
            "sequences": str(tmp_path / "raw" / "sequences"),
            "processed_sequences": str(tmp_path / "processed" / "sequences"),
            "processed_alignments": str(tmp_path / "processed" / "alignments"),
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
            "local_manifest": "",
            "cache_subdir": "uniprot",
            "polyprotein_output": "denv_polyproteins.fasta",
            "e_protein_output": "denv_e_proteins.fasta",
            "manifest_output": "sequence_manifest.tsv",
            "summary_output": "sequence_summary.json",
            "e_protein": {
                "start_1_based": 1,
                "end_1_based_inclusive": 5,
            },
            "validation": {
                "min_e_protein_length": 4,
                "max_e_protein_length": 8,
            },
        },
        "alignment": {
            "method": "python",
            "output_fasta": "denv_e_proteins.aligned.fasta",
            "conservation_output": "conservation.tsv",
            "conservation_by_serotype_output": "conservation_by_serotype.tsv",
            "summary_output": "alignment_summary.json",
            "reference_serotype": "DENV-1",
            "reference_sequence_id": "",
            "python_fallback_max_records": 20,
        },
        "external_tools": {
            "mafft": {"executable": "mafft"},
            "clustal_omega": {"executable": "clustalo"},
        },
    }


class Phase2AlignmentTests(unittest.TestCase):
    def test_external_alignment_commands_are_constructed(self) -> None:
        self.assertEqual(
            build_mafft_command("mafft", "input.fasta"),
            ["mafft", "--auto", "input.fasta"],
        )
        self.assertEqual(
            build_clustal_omega_command("clustalo", "input.fasta", "out.fasta"),
            ["clustalo", "-i", "input.fasta", "-o", "out.fasta", "--force", "--outfmt=fasta"],
        )

    def test_reference_guided_fallback_preserves_insertions(self) -> None:
        aligned_reference, aligned_query = needleman_wunsch("ACGT", "ACGGT")
        self.assertIn("-", aligned_reference)
        self.assertEqual(aligned_reference.replace("-", ""), "ACGT")
        self.assertEqual(aligned_query.replace("-", ""), "ACGGT")

        records = reference_guided_msa(
            [
                FastaRecord("DENV-1|E|0001", "DENV-1|E|0001", "ACGT"),
                FastaRecord("DENV-2|E|0001", "DENV-2|E|0001", "ACGGT"),
                FastaRecord("DENV-3|E|0001", "DENV-3|E|0001", "ACTT"),
            ]
        )
        lengths = {len(record.sequence) for record in records}
        self.assertEqual(len(lengths), 1)
        self.assertTrue(any(record.sequence.count("-") for record in records))
        self.assertIn("-", records[0].sequence)

    def test_conservation_tables_cover_gaps_insertions_and_serotypes(self) -> None:
        records = [
            FastaRecord("DENV-1|E|0001", "DENV-1|E|0001", "A-C"),
            FastaRecord("DENV-2|E|0001", "DENV-2|E|0001", "ATC"),
            FastaRecord("DENV-3|E|0001", "DENV-3|E|0001", "A-C"),
            FastaRecord("DENV-4|E|0001", "DENV-4|E|0001", "AGC"),
        ]
        metadata = {
            "DENV-1|E|0001": AlignmentRecordMetadata(
                "DENV-1|E|0001", "DENV-1", ("A1",), 1, 281, 283
            ),
            "DENV-2|E|0001": AlignmentRecordMetadata(
                "DENV-2|E|0001", "DENV-2", ("A2", "A2b"), 2, 281, 283
            ),
            "DENV-3|E|0001": AlignmentRecordMetadata(
                "DENV-3|E|0001", "DENV-3", ("A3",), 1, 281, 283
            ),
            "DENV-4|E|0001": AlignmentRecordMetadata(
                "DENV-4|E|0001", "DENV-4", ("A4",), 1, 281, 283
            ),
        }

        global_rows, serotype_rows, summary = compute_conservation_tables(
            records,
            metadata,
            reference_id="DENV-1|E|0001",
        )

        self.assertEqual(len(global_rows), 3)
        self.assertEqual(global_rows[1]["reference_e_position"], "")
        self.assertEqual(global_rows[1]["support_count"], 3)
        self.assertEqual(global_rows[1]["gap_fraction"], "0.400000")
        self.assertIn("DENV-2:unique=1", global_rows[1]["serotype_variability"])
        self.assertEqual(summary["reference_columns_numbered"], 2)
        self.assertEqual(summary["insertions_relative_to_reference"], 1)
        self.assertEqual({row["serotype"] for row in serotype_rows}, {"DENV-1", "DENV-2", "DENV-3", "DENV-4"})

    def test_cli_build_alignment_exports_phase2_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "bulk.fasta"
            fasta_path.write_text(
                ">ACC1 DENV-1 envelope isolate\nACGT\n"
                ">ACC2 DENV-2 envelope isolate\nACGGT\n"
                ">ACC3 DENV-3 envelope isolate\nACTT\n"
                ">ACC4 DENV-4 envelope isolate\nACGT\n",
                encoding="utf-8",
            )
            config = phase2_config(tmp_path, fasta_path)
            build_sequence_dataset(config, fetch_remote=False)

            result = run_stage("build-alignment", config, dry_run=False)

            self.assertEqual(result["status"], "completed")
            outputs = result["execution"]["outputs"]
            for path in outputs.values():
                self.assertTrue(Path(path).exists())
            self.assertEqual(result["execution"]["alignment_method"], "python_reference_guided")
            self.assertGreater(result["execution"]["insertions_relative_to_reference"], 0)

            with Path(outputs["conservation"]).open(encoding="utf-8", newline="") as handle:
                conservation_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(any(row["reference_e_position"] == "" for row in conservation_rows))
            self.assertTrue(all(row["support_count"] for row in conservation_rows))


if __name__ == "__main__":
    unittest.main()
