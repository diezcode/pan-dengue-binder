"""Markdown report helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def phase0_report(config: dict[str, Any], stage_results: list[dict[str, Any]]) -> str:
    project = config.get("project", {}).get("name", "pan-dengue-binder")
    lines = [
        "# Phase 0 Status Report",
        "",
        f"Project: `{project}`",
        "",
        "Scope: computational, non-infectious research software foundation.",
        "",
        "## Dry-Run Stages",
        "",
    ]
    for result in stage_results:
        lines.append(f"- `{result['stage']}`: {result['status']}")
    lines.extend(
        [
            "",
            "## Safety Boundary",
            "",
            "This report records software planning and provenance only. It does not validate any medical, experimental, or clinical claim.",
            "",
        ]
    )
    return "\n".join(lines)


def write_report(path: str | Path, content: str) -> Path:
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(content.rstrip() + "\n", encoding="utf-8")
    return output_path

