# Pan-Dengue Binder

Safety-conscious computational research tooling for dengue protein sequence and structure analysis.

This repository is being organized as a reproducible software project for exploring dengue protein conservation, structure mapping, target-atlas reporting, and candidate prioritization hypotheses. It is not a clinical product, medical claim, wet-lab protocol collection, or proof of efficacy.

## Project Scope

The project focuses on safe computational work:

- collecting and documenting dengue sequence and structure metadata;
- comparing DENV-1, DENV-2, DENV-3, and DENV-4 references and broader datasets;
- mapping conservation onto protein structures;
- generating auditable target-atlas reports;
- building dry-run wrappers around external computational tools;
- recording provenance, configuration, and run metadata.

The repository should not contain infectious-virus methods, wet-lab procedures, clinical dosing guidance, or claims that a candidate prevents, treats, or cures disease. Any experimental or translational work would require qualified institutional partners and appropriate review.

## Current Status

Phase 0 foundation is in progress:

- `config/default.yaml` stores project paths, data-source IDs, dry-run defaults, and external tool placeholders.
- `src/pandengue/` contains import-safe Python modules.
- root-level legacy scripts are compatibility wrappers and no longer perform network or external-tool work on import.
- `pandengue` CLI stages can write dry-run metadata.
- `docs/implementation_plan.md` describes the staged implementation roadmap.
- `docs/abbreviations.md` defines common terminology such as ADE, DENV, E protein, EDIII, EDE, and MSA.

Phase 1 sequence data-layer support is also available:

- `pandengue fetch-sequences` can normalize local bulk FASTA files and cached/fetched UniProt references.
- The normalized E-protein dataset is written to `data/processed/sequences/denv_e_proteins.fasta`.
- `data/processed/sequences/sequence_manifest.tsv` records accession, source, serotype, provenance fields, extraction method, duplicate grouping, and filter status.
- `reports/sequence_data.md` states sequence counts per serotype and all filtering outcomes.
- Alignment and target-atlas stages require the manifest before non-dry-run conservation outputs can proceed.

Phase 2 alignment/conservation support is available:

- `pandengue build-alignment` reads the normalized E FASTA and manifest, then writes an aligned FASTA under `data/processed/alignments/`.
- MAFFT or Clustal Omega can be used when configured and installed; a deterministic reference-guided Python fallback supports small test datasets.
- `outputs/target_atlas/conservation.tsv` records alignment-column metrics, entropy, gap fraction, consensus residues, support counts, and reference E/polyprotein numbering.
- `outputs/target_atlas/conservation_by_serotype.tsv` records the same style of metrics per serotype.
- `outputs/target_atlas/alignment_summary.json` records method, parser, reference, counts, and checksums.

Phase 4 target-atlas support is available:

- `pandengue build-target-atlas` ranks epitope hypotheses from Phase 2 conservation data and optional `outputs/target_atlas/structure_mapping.tsv`.
- Candidate patches are filtered by configurable conservation, entropy, gap, and support thresholds, then clustered by sequence proximity and optional 3D residue distance.
- `outputs/target_atlas/epitope_candidates.tsv` and `outputs/target_atlas/ranked_epitopes.json` record transparent score components for conservation, structure, literature/context, ADE context, designability, and template support.
- `reports/target_atlas.md` explains ranked epitopes, failure modes, visualizable residue IDs, and whether 2-4 hypotheses pass the pilot-design decision gate.
- Domain III is treated as one hypothesis among several, not as a hardcoded final target.

## Repository Layout

```text
config/
  default.yaml
docs/
  abbreviations.md
  biomedical_roadmap.md
  implementation_plan.md
src/
  pandengue/
tests/
reports/
```

The current root-level scripts remain for compatibility:

- `fetch_dengue.py`
- `core_alignment_engine.py`
- `find_conserved_target.py`
- `extract_3d_target.py`
- `parse_to_mpnn.py`
- `run_rfdiffusion_wrapper.py`

They now delegate to the package-level planning code.

## CLI Usage

From the repository root, after installing the package in a Python environment:

```powershell
python -m pip install -e .
pandengue --dry-run fetch-sequences
pandengue --dry-run build-alignment
pandengue --dry-run prepare-structure
pandengue --dry-run report
```

In this Codex environment, Python may not be on `PATH`. The bundled runtime can still run the package by setting `PYTHONPATH=src` or using an editable install in a local environment.

Example without installation:

```powershell
$env:PYTHONPATH = "src"
python -m pandengue.cli --dry-run report
```

## Phase Roadmap

The detailed implementation plan is in [docs/implementation_plan.md](docs/implementation_plan.md).

High-level phases:

1. Reproducible repo foundation.
2. Sequence manifest and alignment pipeline.
3. Structure parsing and residue mapping.
4. Target-atlas generation and ranking.
5. Dry-run and then configured external computational design wrappers.
6. Prediction parsing, scoring, clustering, and reporting.
7. Candidate dossiers for qualified review.

## Development Notes

The Phase 0 code is standard-library first so the repository can be tested even before bioinformatics dependencies are installed. Later phases may use optional dependencies such as Biopython, Requests, and PyYAML.

Run unit tests with:

```powershell
python -m unittest discover -s tests
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE).
