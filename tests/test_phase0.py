from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from pandengue.cli import STAGE_PLANNERS, run_report, run_stage
from pandengue.config import load_config
from pandengue.conservation import residue_frequencies, shannon_entropy
from pandengue.metadata import collect_run_metadata
from pandengue.proteinmpnn import build_parse_command
from pandengue.rfdiffusion import build_rfdiffusion_command
from pandengue.sequences import conserved_pattern, find_shared_windows, parse_fasta
from pandengue.structures import count_atom_site_rows


class Phase0Tests(unittest.TestCase):
    def test_default_config_loads_without_pyyaml(self) -> None:
        config = load_config(ROOT / "config" / "default.yaml")
        self.assertEqual(config["project"]["name"], "pan-dengue-binder")
        self.assertTrue(config["pipeline"]["dry_run_default"])
        self.assertEqual(
            config["data_sources"]["uniprot_reference_polyproteins"]["DENV-1"],
            "P33478",
        )

    def test_fasta_parsing_and_conservation(self) -> None:
        records = parse_fasta(
            ">DENV-1 sample\nACDE\n"
            ">DENV-2 sample\nACNE\n"
            ">DENV-3 sample\nACDE\n"
        )
        self.assertEqual([record.identifier for record in records], ["DENV-1", "DENV-2", "DENV-3"])
        self.assertEqual(conserved_pattern(records), "AC-E")

    def test_shared_window_finder(self) -> None:
        windows = find_shared_windows(
            "AAABBBCCCDD",
            ["ZZBBBCCCYY", "XXBBBCCCAA"],
            min_length=3,
        )
        self.assertEqual(windows, ["BBBCCC"])

    def test_conservation_metrics(self) -> None:
        self.assertEqual(shannon_entropy("AAAA"), 0.0)
        frequencies = residue_frequencies("AABC")
        self.assertAlmostEqual(frequencies["A"], 0.5)
        self.assertAlmostEqual(frequencies["B"], 0.25)

    def test_external_commands_are_constructed_not_executed(self) -> None:
        rf_command = build_rfdiffusion_command(
            run_inference_script="",
            target_pdb="target.pdb",
            output_prefix="outputs/design",
            num_designs=3,
            binder_length=55,
        )
        self.assertIn("<RFdiffusion>/scripts/run_inference.py", rf_command)
        self.assertIn("inference.num_designs=3", rf_command)

        mpnn_command = build_parse_command(
            parse_script="",
            input_path="outputs/rfdiffusion",
            output_path="outputs/proteinmpnn/parsed_pdbs.jsonl",
        )
        self.assertIn("<ProteinMPNN>/helper_scripts/parse_multiple_chains.py", mpnn_command)
        self.assertIn("--input_path=outputs/rfdiffusion", mpnn_command)

    def test_structure_counter_handles_mmcif_like_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "mini.cif"
            path.write_text(
                "data_demo\n"
                "ATOM 1 C CA . GLY A 1 1 ? 0 0 0 1.00 1.00 ? 1 GLY A CA 1\n"
                "HETATM 2 O O . HOH A 2 2 ? 0 0 0 1.00 1.00 ? 2 HOH A O 1\n",
                encoding="utf-8",
            )
            self.assertEqual(count_atom_site_rows(path), 2)

    def test_cli_stage_writes_metadata(self) -> None:
        config = load_config(ROOT / "config" / "default.yaml")
        result = run_stage("fetch-sequences", config, dry_run=True)
        self.assertEqual(result["stage"], "fetch-sequences")
        metadata_path = Path(result["metadata_path"])
        self.assertTrue(metadata_path.exists())
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        self.assertEqual(metadata["stage"], "fetch-sequences")
        self.assertTrue(metadata["dry_run"])

    def test_report_command_writes_report(self) -> None:
        config = load_config(ROOT / "config" / "default.yaml")
        result = run_report(config, dry_run=True)
        report_path = Path(result["report_path"])
        self.assertTrue(report_path.exists())
        self.assertIn("Phase 0 Status Report", report_path.read_text(encoding="utf-8"))

    def test_all_planned_stages_are_registered(self) -> None:
        expected = {
            "fetch-sequences",
            "build-alignment",
            "build-target-atlas",
            "prepare-structure",
            "design-backbones",
            "design-sequences",
            "predict-complexes",
            "score-candidates",
        }
        self.assertEqual(set(STAGE_PLANNERS), expected)

    def test_run_metadata_contains_provenance_fields(self) -> None:
        config = load_config(ROOT / "config" / "default.yaml")
        metadata = collect_run_metadata(config, stage="test", dry_run=True)
        self.assertEqual(metadata["stage"], "test")
        self.assertIn("config_hash", metadata)
        self.assertIn("python", metadata)

    def test_legacy_scripts_are_import_safe(self) -> None:
        for script_name in [
            "fetch_dengue.py",
            "core_alignment_engine.py",
            "find_conserved_target.py",
            "extract_3d_target.py",
            "parse_to_mpnn.py",
            "run_rfdiffusion_wrapper.py",
        ]:
            with self.subTest(script=script_name):
                path = ROOT / script_name
                module_name = script_name.replace(".py", "_import_test")
                spec = importlib.util.spec_from_file_location(module_name, path)
                self.assertIsNotNone(spec)
                self.assertIsNotNone(spec.loader)
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)


if __name__ == "__main__":
    unittest.main()

