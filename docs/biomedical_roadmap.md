# Pan-Dengue Binder Biomedical Roadmap

Draft goal: turn this repository into a credible, safety-conscious R&D plan for a pan-serotype dengue intervention. This is not a claim of a cure; it is a staged research plan to discover and validate a candidate that could help solve dengue disease by neutralizing DENV-1, DENV-2, DENV-3, and DENV-4 while minimizing antibody-dependent enhancement (ADE) risk.

## Sources and evidence base to maintain

Keep a living bibliography from:

- WHO dengue fact sheet: https://www.who.int/news-room/fact-sheets/detail/dengue-and-severe-dengue
- CDC dengue clinical information: https://www.cdc.gov/dengue/
- UniProt DENV reference polyproteins, including P33478, P29990, P27915, P09866.
- RCSB PDB dengue E-protein and virion structures, including 1OAN and related antibody/E-protein complexes.
- PubMed literature on dengue broadly neutralizing antibodies, E-dimer epitopes, Domain III, fusion-loop antibodies, quaternary epitopes, and ADE.
- Tool documentation for RFdiffusion, ProteinMPNN, AlphaFold/ColabFold, Rosetta, Biopython, MAFFT/Clustal Omega.

## North-star therapeutic hypothesis

A non-replicating engineered biologic binder can be designed to bind a conserved, surface-accessible, functionally constrained dengue E-protein epitope across all four serotypes. If the binder lacks an antibody Fc function, or uses an Fc-silent/controlled format, it may reduce ADE risk compared with conventional sub-neutralizing antibodies.

Candidate formats to evaluate:

1. De novo miniprotein binder.
2. Nanobody-like binder, preferably Fc-silent if fused.
3. DARPin/affibody/monobody scaffold.
4. Multivalent decoy/trap that crosslinks E proteins or locks the virion in a non-fusogenic state.
5. Later: antibody-derived candidates only with explicit Fc engineering and ADE de-risking.

## Strategic direction

The current repo targets Domain III. That is a reasonable starting point, but the project must not assume Domain III is the only or best target. The next phase should compare multiple target classes:

| Target class | Why interesting | Main risk |
|---|---|---|
| E Domain III lateral ridge / receptor-binding-related surfaces | Known neutralizing antibody target; structurally defined | May be less conserved/accessibility varies |
| E dimer epitope / quaternary epitopes | Often associated with potent broadly neutralizing antibodies | Requires intact virion-like geometry, harder to model from monomer |
| Fusion loop / E stem regions | Functionally constrained and conserved | Some antibodies are cross-reactive but weak/ADE-prone |
| E dimer interface pockets | Could lock prefusion dimer or block rearrangement | Accessibility and strain-state dependence |
| NS1 | Relevant to pathogenesis and diagnostics/therapeutics | Not direct virion neutralization; different disease mechanism |

Recommended primary computational track: **E-protein epitope binder**.
Recommended parallel literature track: **NS1 disease-modifying biologics** as an alternate therapeutic angle.

## Phase 0 — Reframe the repo into a reproducible research project

Deliverables:

- `requirements.txt` or `environment.yml`.
- `config.yaml` with UniProt IDs, PDB IDs, chain IDs, residue ranges, external tool paths, and output folders.
- `data/raw`, `data/processed`, `outputs`, `reports` folders.
- CLI entry points for each stage.
- A report generator that records versions, inputs, parameters, and results.

Immediate code tasks:

1. Make all scripts import-safe using `if __name__ == "__main__"`.
2. Replace hardcoded paths like `/path/to/RFdiffusion` and `/path/to/ProteinMPNN` with config values.
3. Add dry-run mode for external AI tools.
4. Add tests for sequence fetching, alignment, residue mapping, and PDB/mmCIF parsing.

## Phase 1 — Build the target atlas

Goal: identify conserved, surface-exposed dengue target sites with evidence.

Tasks:

1. Download many DENV E-protein sequences, not only four reference polyproteins.
   - Include all four serotypes.
   - Include genotypes and geographically diverse clinical isolates.
2. Perform multiple sequence alignment with MAFFT or Clustal Omega.
3. Compute per-position conservation, entropy, and serotype-specific variability.
4. Map UniProt/polyprotein positions to E-protein positions and PDB/mmCIF residue numbering.
5. Gather structural templates from RCSB:
   - isolated E proteins,
   - E dimers,
   - mature virion structures,
   - antibody-bound structures,
   - serotype-specific structures.
6. Estimate surface exposure/accessibility of conserved residues on biologically relevant structures.
7. Rank epitopes by:
   - conservation,
   - solvent accessibility,
   - functional constraint,
   - antibody/neutralization literature support,
   - ADE risk context,
   - designability for de novo binders.

Decision gate 1: choose 2–4 epitope hypotheses, not one.

## Phase 2 — Design binder candidates computationally

Goal: generate diverse binder designs against prioritized epitopes.

Tasks:

1. Prepare target PDBs with exact chain IDs, residue numbering, protonation/cleanup, and epitope annotations.
2. Generate backbones using RFdiffusion or equivalent.
3. Design sequences using ProteinMPNN.
4. Predict monomer foldability and binder-target complexes using AlphaFold/ColabFold-style methods.
5. Score designs with:
   - predicted interface contacts,
   - buried surface area,
   - shape complementarity,
   - hydrogen bond/salt bridge network,
   - predicted confidence metrics,
   - sequence liability filters,
   - solubility/aggregation risk,
   - immunogenicity risk screens,
   - cross-serotype structural compatibility.
6. Cluster designs to avoid selecting near-duplicates.

Decision gate 2: select ~20–100 computational candidates for non-infectious experimental testing.

## Phase 3 — Non-infectious experimental validation concept

This phase requires qualified wet-lab partners. Keep repo content at planning/specification level.

Recommended assays/concepts:

1. Recombinant E protein or E Domain III binding screen.
2. BLI/SPR affinity and off-rate measurements.
3. Binding to E proteins from all four serotypes.
4. Competition assays against known neutralizing-antibody epitopes.
5. Thermal stability and expression/purification developability.
6. Cell-free or pseudoparticle neutralization models where appropriate.
7. Fc-receptor/ADE de-risking assays only under qualified biosafety and institutional review.

Decision gate 3: advance candidates with broad binding, plausible neutralization proxy activity, and no obvious developability failure.

## Phase 4 — Qualified antiviral and ADE de-risking package

Only certified labs should perform infectious-virus experiments.

Objectives:

1. Confirm neutralization across DENV-1/2/3/4 panels.
2. Test whether sub-neutralizing concentrations enhance infection in relevant Fc receptor contexts.
3. Compare monovalent vs multivalent formats.
4. Evaluate escape risk by analyzing whether target residues are functionally constrained.
5. Begin PK/stability/formulation exploration.

Decision gate 4: choose a lead format and backup series.

## Phase 5 — Translational development sketch

If a lead is validated:

1. Optimize half-life: PEGylation, albumin binding, Fc-silent fusion, or other controlled approaches.
2. Optimize delivery: injectable biologic first; later explore longer-acting formats.
3. Assess immunogenicity and off-target binding.
4. Establish CMC/manufacturing path.
5. Build regulatory pre-IND briefing package.
6. Partner with tropical medicine, virology, and clinical-development groups.

## Key scientific risks

| Risk | Mitigation |
|---|---|
| Target not conserved across real-world isolates | Expand sequence database and update conservation regularly |
| Epitope buried on mature virions | Use mature virion and E-dimer structures, not only monomer structures |
| Designed binder binds recombinant E but not virion | Include quaternary structure models and pseudoparticle/virion-like assays |
| ADE or enhancement-like behavior | Prefer non-Fc binders; test enhancement risk with qualified partners |
| Computational false positives | Use multiple prediction/scoring methods and experimental binding gates |
| Poor manufacturability | Filter early for solubility, stability, aggregation, liabilities |
| Viral escape | Target functionally constrained residues; consider multisite/multivalent designs |

## Immediate next engineering milestones

1. Create `config.yaml` and `requirements.txt`.
2. Replace rough Domain III slicing with real MSA and residue mapping.
3. Add a `target_atlas.py` module that outputs a ranked epitope report.
4. Add mmCIF support for `1OAN.cif` already in the repo.
5. Generate `reports/target_atlas.md` with conservation + structural accessibility evidence.
6. Only then run RFdiffusion design campaigns.

## Ethical and regulatory position

This project should remain computational and non-infectious unless partnered with certified institutions. All disease-treatment claims should be framed as hypotheses until validated by qualified experiments and regulatory review.
