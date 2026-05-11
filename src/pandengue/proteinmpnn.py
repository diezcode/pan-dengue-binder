"""ProteinMPNN command construction."""

from __future__ import annotations


def build_parse_command(
    *,
    parse_script: str,
    input_path: str,
    output_path: str,
) -> list[str]:
    if not parse_script:
        parse_script = "<ProteinMPNN>/helper_scripts/parse_multiple_chains.py"
    return [
        "python",
        parse_script,
        f"--input_path={input_path}",
        f"--output_path={output_path}",
    ]


def plan_proteinmpnn(config: dict) -> dict:
    tool = config.get("external_tools", {}).get("proteinmpnn", {})
    command = build_parse_command(
        parse_script=str(tool.get("parse_multiple_chains_script") or ""),
        input_path="outputs/rfdiffusion/",
        output_path="outputs/proteinmpnn/parsed_pdbs.jsonl",
    )
    return {
        "action": "construct_proteinmpnn_parse_command",
        "dry_run_safe": True,
        "command": command,
        "outputs": ["outputs/proteinmpnn/"],
    }

