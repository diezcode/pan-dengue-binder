"""RFdiffusion command construction.

Phase 0 only builds auditable dry-run commands. Later phases can add execution
once external tools are installed and configured.
"""

from __future__ import annotations

from pathlib import Path


def build_rfdiffusion_command(
    *,
    run_inference_script: str,
    target_pdb: str,
    output_prefix: str,
    num_designs: int,
    binder_length: int,
    target_chain: str = "A",
    target_start: int = 1,
    target_end: int = 100,
) -> list[str]:
    if not run_inference_script:
        run_inference_script = "<RFdiffusion>/scripts/run_inference.py"
    contig = f"[{target_chain}{target_start}-{target_end}/{binder_length}-{binder_length}]"
    return [
        "python",
        run_inference_script,
        f"inference.output_prefix={output_prefix}",
        f"inference.input_pdb={target_pdb}",
        f"inference.num_designs={num_designs}",
        f"contigmap.contigs={contig}",
    ]


def plan_rfdiffusion(config: dict) -> dict:
    defaults = config.get("design_defaults", {})
    tool = config.get("external_tools", {}).get("rfdiffusion", {})
    command = build_rfdiffusion_command(
        run_inference_script=str(tool.get("run_inference_script") or ""),
        target_pdb="data/processed/structures/<target>.pdb",
        output_prefix=str(defaults.get("output_prefix", "outputs/rfdiffusion/design")),
        num_designs=int(defaults.get("num_designs", 10)),
        binder_length=int(defaults.get("binder_length", 65)),
        target_chain=str(defaults.get("target_chain", "A")),
    )
    return {
        "action": "construct_rfdiffusion_design_command",
        "dry_run_safe": True,
        "command": command,
        "outputs": ["outputs/rfdiffusion/"],
    }


def command_as_text(command: list[str]) -> str:
    return " ".join(str(part) for part in command)


def validate_existing_file(path: str) -> Path:
    resolved = Path(path)
    if not resolved.exists():
        raise FileNotFoundError(path)
    return resolved

