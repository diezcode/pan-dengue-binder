"""Small conservation metrics used by Phase 0 tests and later atlas work."""

from __future__ import annotations

import math
from collections import Counter


def shannon_entropy(residues: str) -> float:
    if not residues:
        return 0.0
    counts = Counter(residues)
    total = len(residues)
    return -sum((count / total) * math.log2(count / total) for count in counts.values())


def residue_frequencies(residues: str) -> dict[str, float]:
    if not residues:
        return {}
    counts = Counter(residues)
    total = len(residues)
    return {residue: count / total for residue, count in sorted(counts.items())}

