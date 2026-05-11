"""Run metadata and provenance helpers."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import repo_root


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, separators=(",", ":"))
    return sha256_text(payload)


def _git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root(),
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except Exception:
        return None
    commit = result.stdout.strip()
    return commit or None


def collect_run_metadata(
    config: dict[str, Any],
    *,
    stage: str,
    dry_run: bool,
    inputs: dict[str, str] | None = None,
    tools: dict[str, str] | None = None,
) -> dict[str, Any]:
    return {
        "stage": stage,
        "dry_run": dry_run,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": config.get("project", {}).get("name", "pan-dengue-binder"),
        "config_hash": config_hash(config),
        "git_commit": _git_commit(),
        "python": {
            "version": sys.version.split()[0],
            "implementation": platform.python_implementation(),
        },
        "platform": platform.platform(),
        "inputs": inputs or {},
        "tools": tools or {},
    }


def write_run_metadata(path: str | Path, metadata: dict[str, Any]) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return output_path

