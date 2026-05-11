"""Compatibility wrapper for ProteinMPNN input planning."""

from __future__ import annotations

import json
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from pandengue.cli import run_stage
from pandengue.config import load_config
from pandengue.proteinmpnn import build_parse_command


def prepare_for_mpnn(rfdiffusion_output_dir: str, mpnn_input_dir: str) -> list[str]:
    """Return the planned parser command without executing external tools."""
    return build_parse_command(
        parse_script="",
        input_path=rfdiffusion_output_dir,
        output_path=str(Path(mpnn_input_dir) / "parsed_pdbs.jsonl"),
    )


def main() -> int:
    result = run_stage("design-sequences", load_config(), dry_run=True)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

