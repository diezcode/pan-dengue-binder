"""Compatibility wrapper for the sequence-fetch planning stage."""

from __future__ import annotations

import json
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from pandengue.cli import run_stage
from pandengue.config import load_config


def main() -> int:
    result = run_stage("fetch-sequences", load_config(), dry_run=True)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

