"""Compatibility wrapper for RFdiffusion command planning."""

from __future__ import annotations

import json
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from pandengue.cli import run_stage
from pandengue.config import load_config
from pandengue.rfdiffusion import build_rfdiffusion_command


def run_rfdiffusion(
    target_pdb: str,
    output_prefix: str,
    num_designs: int,
    binder_length: int,
) -> list[str]:
    """Return the planned RFdiffusion command without executing it."""
    return build_rfdiffusion_command(
        run_inference_script="",
        target_pdb=target_pdb,
        output_prefix=output_prefix,
        num_designs=num_designs,
        binder_length=binder_length,
    )


def main() -> int:
    result = run_stage("design-backbones", load_config(), dry_run=True)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

