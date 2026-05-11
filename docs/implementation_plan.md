# Pan-Dengue Binder Implementation Plan

Goal: turn this repository into a reproducible, safety-conscious computational pipeline for discovering a pan-serotype dengue biologic candidate. The intended product is not a claimed cure. The intended product is a ranked set of non-infectious, experimentally testable binder candidates that may neutralize DENV-1, DENV-2, DENV-3, and DENV-4 while minimizing antibody-dependent enhancement (ADE) risk.

## North Star

Build a computational discovery engine for a non-replicating biologic binder that:

- binds conserved, surface-accessible, functionally constrained dengue E-protein epitopes across DENV-1/2/3/4;
- avoids Fc receptor mediated ADE risk by defaulting to miniprotein, DARPin, monobody, affibody, or nanobody-like formats without an active Fc domain;
- resists viral escape by targeting one or more structurally constrained sites rather than relying on a single reference-sequence motif;
- produces auditable outputs that can be handed to qualified wet-lab partners for non-infectious binding and neutralization-proxy testing.

## Safety and Scope

This repository should remain computational and non-infectious. It may define experimental acceptance criteria, but it should not provide instructions for culturing, propagating, concentrating, engineering, or infecting with dengue virus. Infectious-virus neutralization and ADE de-risking must be performed only by qualified institutional partners under appropriate biosafety, ethics, and regulatory oversight.

Every user-facing claim should use research language: "candidate", "hypothesis", "predicted", "screened", "requires validation". Avoid "cure", "guaranteed", "100% immutable", or "virus cannot adapt".

## Current Implementation Status

As of the current repository state:

| Phase | Status | Notes |
|---|---|---|
| Phase 0: Reframe and Reproducibility | Implemented | Package layout, config, import-safe wrappers, CLI planning/execution hooks, dry-run metadata, and first tests are present. |
| Phase 1: Sequence Data Layer | Implemented | Manifest-backed sequence ingestion, E-protein normalization, deduplication, dataset gates, and sequence report outputs are present. |
| Phase 2: Alignment and Conservation Engine | Implemented | Alignment execution/fallback, conservation metrics, serotype metrics, reference numbering, and target-atlas conservation outputs are present. |
| Phase 3: Structure and Residue Mapping | Not implemented yet | Full residue mapping, structure-aware accessibility, glycan/context annotations, and `outputs/target_atlas/structure_mapping.tsv` generation still need to be built. |
| Phase 4: Target Atlas and Epitope Ranking | Implemented | Candidate patch generation, clustering, scoring, ranked outputs, decision gate, and `reports/target_atlas.md` are present. Uses optional Phase 3 structure mapping when available. |
| Phase 5 and later | Not implemented yet | Binder design, sequence design, prediction/scoring, ADE package, and final dossier/report gates remain future work. |

## Target Architecture

```text
config/
  default.yaml
data/
  raw/
    sequences/
    structures/
  processed/
    alignments/
    mappings/
    structures/
outputs/
  target_atlas/
  rfdiffusion/
  proteinmpnn/
  predictions/
  scored_candidates/
reports/
  target_atlas.md
  design_campaign.md
src/pandengue/
  cli.py
  config.py
  logging.py
  sequences.py
  alignments.py
  conservation.py
  structures.py
  residue_mapping.py
  epitope_atlas.py
  rfdiffusion.py
  proteinmpnn.py
  colabfold.py
  scoring.py
  reporting.py
tests/
```

The existing root-level scripts should become thin wrappers or be migrated into `src/pandengue/`. The pipeline should be runnable stage-by-stage and end-to-end.

## Phase 0: Reframe and Reproducibility

Implementation status: implemented.

Deliverables:

- `README.md` rewritten to describe a computational candidate-discovery project, not a cure.
- `pyproject.toml` or `requirements.txt` with pinned core dependencies.
- `config/default.yaml` containing data sources, serotype labels, PDB IDs, chain IDs, external tool paths, output folders, scoring thresholds, and dry-run settings.
- import-safe Python modules with `if __name__ == "__main__"` only in CLI entry points.
- deterministic run metadata: timestamp, git commit if available, config hash, input file checksums, tool versions.

Implementation tasks:

1. Move logic from `fetch_dengue.py`, `find_conserved_target.py`, `core_alignment_engine.py`, `extract_3d_target.py`, `parse_to_mpnn.py`, and `run_rfdiffusion_wrapper.py` into modules.
2. Add a CLI:
   - `pandengue fetch-sequences`
   - `pandengue build-alignment`
   - `pandengue build-target-atlas`
   - `pandengue prepare-structure`
   - `pandengue design-backbones`
   - `pandengue design-sequences`
   - `pandengue predict-complexes`
   - `pandengue score-candidates`
   - `pandengue report`
3. Add `--dry-run` to every stage that invokes external tools.
4. Add tests for config loading, FASTA parsing, simple conservation scoring, mmCIF parsing, residue mapping, and command construction.

Done when:

- all existing scripts can be imported without side effects;
- a dry-run end-to-end command produces a run folder and report without requiring RFdiffusion or ProteinMPNN;
- tests pass locally.

## Phase 1: Sequence Data Layer

Implementation status: implemented.

Goal: replace four-reference-sequence assumptions with a broad, auditable sequence set.

Data sources:

- UniProt dengue reference polyproteins for stable reference anchors.
- NCBI Virus or NCBI Datasets for broader geographically and temporally diverse DENV sequences.
- A local manifest recording accession, serotype, host, country/region, collection date if available, source database, and retrieval date.

Implementation tasks:

1. Add `sequences.py` to fetch/cache FASTA and metadata.
2. Support both curated references and bulk isolate data.
3. Normalize records into:
   - `data/processed/sequences/denv_e_proteins.fasta`
   - `data/processed/sequences/sequence_manifest.tsv`
4. Extract E protein and optionally NS1 using annotation-aware coordinates where available. Fall back to reference-guided mapping only when annotations are missing.
5. De-duplicate identical sequences but keep metadata counts.
6. Enforce minimum dataset gates:
   - all four serotypes represented;
   - configurable minimum number of E sequences per serotype;
   - all sequence filters written to the report.

Done when:

- the report states exactly how many sequences per serotype were used;
- every sequence has an accession and source;
- no conservation result can be produced without a manifest.

## Phase 2: Alignment and Conservation Engine

Implementation status: implemented.

Goal: identify conserved regions with real multiple sequence alignment and entropy metrics.

Implementation tasks:

1. Add MAFFT or Clustal Omega integration with a pure-Python fallback for small tests.
2. Parse alignments with Biopython `AlignIO`.
3. Compute per-position:
   - residue frequency by serotype and globally;
   - Shannon entropy;
   - gap fraction;
   - consensus residue;
   - serotype-specific variability;
   - real-world isolate support count.
4. Map alignment columns back to reference E-protein numbering and polyprotein numbering.
5. Export:
   - `outputs/target_atlas/conservation.tsv`
   - `outputs/target_atlas/conservation_by_serotype.tsv`
   - `outputs/target_atlas/alignment_summary.json`

Done when:

- conserved positions are derived from alignment columns, not string containment;
- each conserved residue can be traced back to reference numbering and supporting isolates;
- tests cover gaps, insertions, and mixed serotype records.

## Phase 3: Structure and Residue Mapping

Implementation status: not implemented yet.

Goal: stop using rough residue slices and map sequence conservation onto real dengue structures.

Implementation tasks:

1. Add mmCIF-first structure parsing for PDB entries such as `1OAN.cif` and related E-protein, E-dimer, virion, and antibody-bound structures.
2. Use RCSB metadata to capture entity IDs, chain IDs, experimental method, resolution where relevant, organism/serotype, ligand/glycan annotations, and biological assembly.
3. Build `residue_mapping.py`:
   - UniProt/polyprotein position to mature E-protein position;
   - E-protein position to PDB/mmCIF auth and label residue IDs;
   - chain-aware mappings;
   - insertion code handling.
4. Compute or import structural features:
   - solvent accessible surface area;
   - distance to glycans;
   - proximity to known antibody epitopes;
   - whether the residue is exposed in monomer, dimer, and virion-like context;
   - whether the site is likely buried or artificially exposed by cropping.
5. Export cleaned target structures with annotated epitope residues:
   - `data/processed/structures/<target_id>.pdb`
   - `outputs/target_atlas/structure_mapping.tsv`

Done when:

- `extract_3d_target.py` is replaced by structure-aware target preparation;
- every design target records exact chain and residue IDs;
- the report distinguishes recombinant-domain accessibility from mature-virion accessibility.

## Phase 4: Target Atlas and Epitope Ranking

Implementation status: implemented. This phase can consume optional Phase 3 structure mapping if `outputs/target_atlas/structure_mapping.tsv` is available, but Phase 3 itself is still open.

Goal: choose 2-4 target hypotheses before running expensive design campaigns.

Target classes to compare:

- E Domain III receptor-binding/lateral-ridge regions.
- E-dimer epitope or EDE-like quaternary surfaces.
- Fusion-loop and surrounding E Domain II surfaces.
- E dimer interface pockets that might lock prefusion E.
- E stem/proximal regions if accessible in relevant structures.
- NS1 only as a separate disease-modifying track, not a direct virion-neutralization track.

Ranking criteria:

- pan-serotype conservation;
- low entropy across real isolates;
- surface accessibility in relevant structures;
- functional constraint and escape cost;
- known neutralization or protective-antibody literature support;
- ADE risk context;
- designability for non-Fc binders;
- glycan and conformational-breathing risk;
- availability of high-quality structural templates.

Implementation tasks:

1. Add `epitope_atlas.py` to generate candidate patches from conserved exposed residues.
2. Cluster nearby residues into structural epitopes.
3. Score and rank patches with transparent weights from config.
4. Produce `reports/target_atlas.md` with:
   - ranked epitopes;
   - rationale;
   - visualizable residue lists;
   - failure modes;
   - recommendation for pilot design.

Decision gate:

- Advance 2-4 epitopes only if they have clear sequence, structural, and literature support.
- If only one epitope passes, do not claim pan-dengue robustness; expand data or structures first.

Done when:

- the repo can explain why a target was chosen;
- Domain III is treated as one hypothesis, not as a hardcoded truth.

## Phase 5: Binder Backbone Design

Goal: generate diverse binder backbones against prioritized epitopes.

Implementation tasks:

1. Replace hardcoded RFdiffusion paths with config.
2. Generate RFdiffusion commands from structured target definitions.
3. Support hotspot residues, cropped target structures, binder length ranges, number of designs, seeds, and noise-scale parameters.
4. Validate that target residues in the contig string exist in the prepared PDB.
5. Save every design campaign as:
   - input config;
   - command lines;
   - target PDB;
   - hotspot residue list;
   - RFdiffusion outputs;
   - run metadata.

Pilot campaign:

- Run 50-200 backbones per epitope in dry-run or small real campaigns.
- Inspect whether binders contact the intended epitope.

Scale campaign:

- Only after the pilot passes, generate 1,000-10,000 backbones per selected target, depending on compute.

Done when:

- no RFdiffusion command can be built from a vague residue range;
- every generated backbone is linked to a specific epitope and config;
- target-chain and binder-chain conventions are consistent.

## Phase 6: Sequence Design

Goal: convert backbones to amino acid sequences while keeping the dengue target fixed.

Implementation tasks:

1. Add `proteinmpnn.py` wrapper around ProteinMPNN helper scripts:
   - parse PDBs to JSONL;
   - assign fixed target chains;
   - optionally make fixed-position dictionaries;
   - run `protein_mpnn_run.py`.
2. Support soluble model settings and multiple sampling temperatures.
3. Add sequence liability filters:
   - extreme charge;
   - long hydrophobic patches;
   - unpaired cysteines unless intentionally designed;
   - N-linked glycosylation motifs if problematic for expression format;
   - protease-sensitive or aggregation-prone motifs where predictable.
4. Export:
   - designed FASTA;
   - binder-target PDB paths;
   - per-sequence ProteinMPNN scores;
   - liability annotations.

Done when:

- designed chains and fixed chains are explicit;
- the target sequence is not redesigned by mistake;
- sequence filters run before expensive structure prediction.

## Phase 7: In-Silico Prediction and Scoring

Goal: reduce computational false positives before any experimental partner work.

Implementation tasks:

1. Add ColabFold/AlphaFold-compatible input generation for binder-target complexes.
2. Parse predicted outputs and score:
   - binder monomer confidence;
   - complex confidence;
   - predicted aligned error across the interface;
   - interface contact count;
   - buried surface area;
   - shape complementarity if available;
   - hydrogen bond and salt-bridge plausibility;
   - predicted off-target/self-binding red flags;
   - cross-serotype structural compatibility.
3. Add optional Rosetta InterfaceAnalyzer scoring for top candidates.
4. Cluster candidates by sequence and structural similarity to avoid near-duplicate picks.
5. Export:
   - `outputs/scored_candidates/candidates.tsv`
   - `outputs/scored_candidates/top_candidates.fasta`
   - `reports/design_campaign.md`

Initial computational pass thresholds should be conservative and config-driven. Example fields:

- `pae_interaction_max`
- `iptm_min`
- `interface_contacts_min`
- `binder_plddt_min`
- `liability_score_max`
- `serotype_coverage_min`

Done when:

- candidates are ranked by multiple independent signals;
- no single model score is treated as proof of efficacy;
- every final candidate has traceable input data and generated artifacts.

## Phase 8: ADE and Translational Risk Package

Goal: make ADE avoidance a first-class design requirement.

Implementation tasks:

1. Add candidate format labels:
   - miniprotein binder;
   - nanobody-like binder;
   - DARPin/monobody/affibody;
   - Fc-fusion only if Fc-silent and explicitly justified.
2. Add ADE risk annotations:
   - Fc present or absent;
   - Fc receptor binding expected or disabled;
   - epitope class associated with neutralization versus weak cross-reactivity;
   - predicted sub-neutralizing binding concern;
   - multivalent avidity risk/benefit note.
3. Add a non-infectious validation brief template for qualified partners:
   - expression and purification feasibility;
   - recombinant E/E-dimer binding;
   - all-serotype binding panel;
   - BLI/SPR affinity and off-rate;
   - competition with known neutralizing epitope probes;
   - pseudoparticle or other qualified neutralization-proxy assay where appropriate.

Done when:

- every candidate has an ADE risk note;
- Fc-active formats are excluded by default;
- the repo clearly says infectious-virus testing is out of scope.

## Phase 9: Reporting, Review, and Release Gates

Goal: make the project auditable enough for external scientific review.

Reports:

- `target_atlas.md`: why these epitopes were chosen.
- `design_campaign.md`: how candidates were generated.
- `candidate_dossier_<id>.md`: one-page dossier per top candidate.
- `safety_scope.md`: what the repo does and does not do.

Each candidate dossier should include:

- target epitope and residue IDs;
- serotype coverage evidence;
- structure template;
- design method and parameters;
- sequence;
- predicted structure metrics;
- liability filters;
- ADE risk position;
- recommended non-infectious validation next step;
- known reasons the candidate might fail.

Release gate for "computational candidate":

- broad target conservation supported by real isolate data;
- target exposed in relevant structures;
- predicted binder-target interface passes thresholds;
- candidate passes sequence/developability filters;
- at least two backup candidates from a different cluster;
- report generated from reproducible run metadata.

## Suggested Implementation Order

1. Build repo foundation: config, package layout, CLI, dry-run, tests.
2. Replace current sequence fetching and conservation scan.
3. Add alignment and conservation report.
4. Add mmCIF parsing and residue mapping.
5. Build target atlas and ranking report.
6. Update RFdiffusion wrapper to use atlas targets.
7. Update ProteinMPNN wrapper to preserve target chains.
8. Add prediction/scoring parsers.
9. Add candidate dossiers and safety report.
10. Only then run a real design campaign.

## Near-Term Sprint Backlog

Sprint 1:

- Create `pyproject.toml`, `config/default.yaml`, and `src/pandengue/`.
- Convert existing scripts into import-safe modules.
- Add dry-run CLI.
- Add first tests.

Sprint 2:

- Implement sequence manifest and cached FASTA retrieval.
- Add broad sequence ingestion from NCBI datasets exports.
- Add alignment runner and conservation metrics.
- Generate first `target_atlas` conservation tables.

Sprint 3:

- Implement mmCIF parsing for `1OAN.cif`.
- Add residue mapping between reference E positions and structure residues.
- Add surface/accessibility annotations.
- Produce the first ranked target atlas report.

Sprint 4:

- Refactor RFdiffusion wrapper around config-driven target definitions.
- Refactor ProteinMPNN bridge with fixed target chain handling.
- Add command construction tests.
- Run dry-run design campaign from atlas output.

Sprint 5:

- Add ColabFold/AlphaFold output parser.
- Add scoring, clustering, and candidate dossier generation.
- Produce a small pilot candidate report.

## Required External Resources

Scientific and clinical context:

- WHO dengue overview: https://www.who.int/mega-menu/health-topics/popular/dengue
- CDC dengue clinical care: https://www.cdc.gov/dengue/hcp/clinical-care/index.html
- WHO dengue vaccine Q&A: https://www.who.int/news-room/questions-and-answers/item/dengue-vaccines
- CDC dengue vaccine information: https://www.cdc.gov/dengue/hcp/vaccine/index.html

Data sources:

- UniProt/EBI proteins API: https://www.ebi.ac.uk/proteins/api
- NCBI Datasets virus downloads: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/virus/virus-download/
- NCBI Virus help: https://www.ncbi.nlm.nih.gov/labs/virus/vssi/docs/help/
- RCSB PDB Data API: https://data.rcsb.org/
- RCSB PDB Search API: https://search.rcsb.org/

Computational tools:

- Biopython AlignIO: https://biopython.org/wiki/AlignIO
- RFdiffusion: https://github.com/RosettaCommons/RFdiffusion
- ProteinMPNN: https://github.com/dauparas/ProteinMPNN
- ColabFold: https://github.com/sokrypton/ColabFold
- Rosetta InterfaceAnalyzer: https://docs.rosettacommons.org/docs/latest/application_documentation/analysis/interface-analyzer

Literature anchors:

- ADE and Fc receptor risk: https://www.nature.com/articles/s41577-020-00410-0
- Dengue neutralizing antibody targets and ADE review: https://pmc.ncbi.nlm.nih.gov/articles/PMC10272415/
- Dengue monoclonal antibody development review: https://pubmed.ncbi.nlm.nih.gov/40916325/
- E-dimer epitope broadly protective antibody preprint: https://www.medrxiv.org/content/10.1101/2024.04.30.24306574v1
