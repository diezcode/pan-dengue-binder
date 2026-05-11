"""Configuration loading for the computational pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def repo_root() -> Path:
    """Return the repository root for an editable source checkout."""
    return Path(__file__).resolve().parents[2]


def default_config_path() -> Path:
    return repo_root() / "config" / "default.yaml"


def _parse_scalar(value: str) -> Any:
    value = value.strip()
    if value == "":
        return ""
    if value.lower() == "true":
        return True
    if value.lower() == "false":
        return False
    if value.lower() in {"null", "none"}:
        return None
    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def _parse_simple_yaml(text: str) -> dict[str, Any]:
    """Parse the small YAML subset used by config/default.yaml.

    PyYAML is optional for later phases. Phase 0 stays runnable with only the
    standard library, so this parser supports nested mappings and scalar values.
    """
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        if "\t" in line:
            raise ValueError(f"Tabs are not supported in YAML at line {line_number}")

        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()
        if stripped.startswith("- "):
            raise ValueError(
                f"Lists are not supported by the Phase 0 YAML parser at line {line_number}"
            )
        if ":" not in stripped:
            raise ValueError(f"Expected key/value YAML mapping at line {line_number}")

        key, raw_value = stripped.split(":", 1)
        key = key.strip()
        raw_value = raw_value.strip()
        if not key:
            raise ValueError(f"Empty YAML key at line {line_number}")

        while stack and indent <= stack[-1][0]:
            stack.pop()
        if not stack:
            raise ValueError(f"Invalid indentation at line {line_number}")

        parent = stack[-1][1]
        if raw_value == "":
            child: dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = _parse_scalar(raw_value)

    return root


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load project configuration.

    If PyYAML is installed, it is used. Otherwise the repository's simple YAML
    subset is parsed with a standard-library fallback.
    """
    config_path = Path(path) if path is not None else default_config_path()
    text = config_path.read_text(encoding="utf-8")

    try:
        import yaml  # type: ignore
    except ModuleNotFoundError:
        return _parse_simple_yaml(text)

    data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError(f"Config must be a mapping: {config_path}")
    return data


def get_nested(config: dict[str, Any], dotted_key: str, default: Any = None) -> Any:
    current: Any = config
    for part in dotted_key.split("."):
        if not isinstance(current, dict) or part not in current:
            return default
        current = current[part]
    return current


def resolve_path(config: dict[str, Any], dotted_key: str) -> Path:
    raw = get_nested(config, dotted_key)
    if raw is None:
        raise KeyError(f"Missing path setting: {dotted_key}")
    path = Path(str(raw))
    if path.is_absolute():
        return path
    return repo_root() / path

