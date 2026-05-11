"""Command-line interface for the Phase 0 computational pipeline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

from .alignments import build_alignment_dataset, plan_alignment
from .config import get_nested, load_config, resolve_path
from .epitope_atlas import build_target_atlas, plan_target_atlas
from .metadata import collect_run_metadata, write_run_metadata
from .proteinmpnn import plan_proteinmpnn
from .reporting import phase0_report, write_report
from .rfdiffusion import plan_rfdiffusion
from .sequences import build_sequence_dataset, plan_sequence_fetch, require_sequence_manifest
from .structures import plan_structure_preparation


StagePlanner = Callable[[dict[str, Any]], dict[str, Any]]

STAGE_PLANNERS: dict[str, StagePlanner] = {
    "fetch-sequences": plan_sequence_fetch,
    "build-alignment": plan_alignment,
    "build-target-atlas": plan_target_atlas,
    "prepare-structure": plan_structure_preparation,
    "design-backbones": plan_rfdiffusion,
    "design-sequences": plan_proteinmpnn,
    "predict-complexes": lambda config: {
        "action": "prepare_prediction_inputs_and_parse_prediction_outputs",
        "dry_run_safe": True,
        "outputs": ["outputs/predictions/"],
    },
    "score-candidates": lambda config: {
        "action": "score_and_rank_candidates_from_existing_outputs",
        "dry_run_safe": True,
        "outputs": ["outputs/scored_candidates/"],
    },
}

MANIFEST_REQUIRED_STAGES = {"build-alignment", "build-target-atlas"}


def _metadata_path(config: dict[str, Any], stage: str) -> Path:
    outputs_root = resolve_path(config, "paths.outputs")
    filename = str(get_nested(config, "pipeline.run_metadata_filename", "run_metadata.json"))
    return outputs_root / "phase0" / stage / filename


def run_stage(stage: str, config: dict[str, Any], *, dry_run: bool) -> dict[str, Any]:
    if stage not in STAGE_PLANNERS:
        raise KeyError(f"Unknown stage: {stage}")
    plan = STAGE_PLANNERS[stage](config)
    execution: dict[str, Any] | None = None
    status = "planned" if dry_run else "ready_for_implementation"

    if not dry_run:
        if stage == "fetch-sequences":
            execution = build_sequence_dataset(config, fetch_remote=True)
            status = "completed"
        elif stage == "build-alignment":
            execution = build_alignment_dataset(config)
            status = "completed"
        elif stage == "build-target-atlas":
            execution = build_target_atlas(config)
            status = "completed"
        elif stage in MANIFEST_REQUIRED_STAGES:
            manifest_path = require_sequence_manifest(config)
            plan["manifest_path"] = str(manifest_path)

    metadata = collect_run_metadata(
        config,
        stage=stage,
        dry_run=dry_run,
        tools={stage: plan.get("action", "planned")},
    )
    metadata_path = write_run_metadata(_metadata_path(config, stage), metadata)
    result = {
        "stage": stage,
        "status": status,
        "dry_run": dry_run,
        "plan": plan,
        "metadata_path": str(metadata_path),
    }
    if execution is not None:
        result["execution"] = execution
    return result


def run_report(config: dict[str, Any], *, dry_run: bool) -> dict[str, Any]:
    results = [run_stage(stage, config, dry_run=dry_run) for stage in STAGE_PLANNERS]
    report_path = resolve_path(config, "paths.reports") / "phase0_status.md"
    write_report(report_path, phase0_report(config, results))
    return {
        "stage": "report",
        "status": "written",
        "dry_run": dry_run,
        "report_path": str(report_path),
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pandengue",
        description="Safe computational research pipeline utilities.",
    )
    parser.add_argument("--config", default=None, help="Path to a YAML config file.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=None,
        help="Plan the action and write metadata without invoking external tools.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)
    for stage in STAGE_PLANNERS:
        subparsers.add_parser(stage, help=f"Plan or run stage: {stage}")
    subparsers.add_parser("report", help="Write a Phase 0 status report.")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    config = load_config(args.config)
    dry_run = (
        bool(args.dry_run)
        if args.dry_run is not None
        else bool(get_nested(config, "pipeline.dry_run_default", True))
    )

    if args.command == "report":
        result = run_report(config, dry_run=dry_run)
    else:
        result = run_stage(args.command, config, dry_run=dry_run)

    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
