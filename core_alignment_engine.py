"""Compatibility wrapper for the alignment planning stage."""

from __future__ import annotations

import json
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from pandengue.cli import run_stage
from pandengue.config import load_config
from pandengue.sequences import find_shared_windows


def find_immutable_targets(seqs: dict[str, str], min_length: int = 10) -> list[str]:
    """Backward-compatible helper for older notebooks/scripts."""
    reference = seqs["DENV-1"]
    comparisons = [sequence for name, sequence in seqs.items() if name != "DENV-1"]
    return find_shared_windows(reference, comparisons, min_length=min_length)


def main() -> int:
    result = run_stage("build-alignment", load_config(), dry_run=True)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

