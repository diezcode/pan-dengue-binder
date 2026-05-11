"""Structure-file helpers for the Phase 0 foundation."""

from __future__ import annotations

from pathlib import Path


def count_atom_site_rows(path: str | Path) -> int:
    """Count ATOM/HETATM rows in PDB-like or mmCIF-like files."""
    count = 0
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.lstrip()
            if stripped.startswith("ATOM ") or stripped.startswith("HETATM "):
                count += 1
    return count


def plan_structure_preparation(config: dict) -> dict:
    structures = config.get("data_sources", {}).get("rcsb_structures", {})
    return {
        "action": "prepare_structure_metadata_and_mappings",
        "dry_run_safe": True,
        "seed_structure": structures.get("seed_structure"),
        "outputs": [
            "data/processed/structures/",
            "outputs/target_atlas/structure_mapping.tsv",
        ],
    }

