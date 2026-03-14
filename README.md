# Pan-Dengue Binder: Generative AI In-Silico Pipeline

<p align="center">
  <img src="https://img.shields.io/badge/Status-Active%20Development-green.svg" alt="Project Status">
  <img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License">
  <img src="https://img.shields.io/badge/Tech-Python%20%7C%20CUDA%20%7C%20WSL2-yellow.svg" alt="Tech Stack">
  <img src="https://img.shields.io/badge/AI-RFdiffusion%20%7C%20ProteinMPNN%20%7C%20AlphaFold2-purple.svg" alt="AI Models">
</p>

## 🧬 Project Overview & The Context
The goal of this project is to build a modern **in-silico Bio-AI Pipeline** designed for the research and generative development of a universal neutralizer against the Dengue virus. 

Instead of relying on the error-prone human immune system, we approach virology as software engineers. We design a software architecture that orchestrates cutting-edge generative AI models to mathematically design artificial "mini-proteins" (de novo binders) from scratch. These computational molecules are designed to mechanically lock and seal the virus, preventing it from binding to human cells.

### The Biological Problem: The ADE Effect
Dengue fever is notoriously difficult to vaccinate against because the virus exists in four distinct serotypes (variants). The primary challenge in combating Dengue is the **Antibody-Dependent Enhancement (ADE) effect**. 

If a patient is infected with Dengue Type 1, they develop antibodies against it. However, if they are later infected with Type 2, the immune system deploys the old (Type 1) antibodies. These antibodies only partially fit the new variant. Instead of neutralizing the virus, the antibody accidentally acts as a "Trojan Horse", helping the live virus to heavily invade human immune cells (macrophages). This explosive viral replication causes severe, often fatal, hemorrhagic fever.

### The Solution Strategy: The Achilles Heel
To bypass the deadly ADE effect, our pipeline targets a single, highly conserved region on the virus's surface—specifically, **Domain III of the E-Protein**. 

Algorithmic alignment of all four Dengue variants proves that this tiny physical structure is **100% identical across all serotypes**. The virus cannot mutate this spot without destroying its own ability to infect cells. 
Our strategy leverages generative AI to craft an entirely new, hyper-specific molecule (a "binder"). This artificial molecule is engineered to fit perfectly and exclusively onto this immutable weak spot—mechanically sealing the virus without triggering the ADE mechanism.

---

## 🏗️ AI Architecture & Pipeline

This project is not a single script, but a chained assembly line of multiple state-of-the-art neural networks running locally on consumer hardware (e.g., RTX 4080 via WSL2).

### 1. Data Acquisition & Target Isolation (Python Layer)
* **`fetch_dengue.py`**: Interacts with the UniProt biological database to download the raw amino acid sequences of the virus.
* **`core_alignment_engine.py` & `find_conserved_target.py`**: Algorithmically scans the four viral variants, performing Multiple Sequence Alignment (MSA) to map the exact conserved target zones.
* **`extract_3d_target.py`**: Parses the structural 3D CAD files (.pdb) from the Protein Data Bank and slices out only the precise 3D coordinates (X, Y, Z) of the target motif (Domain III). This becomes the invariant template for the AI.

### 2. Generative Design (RFdiffusion)
* **The Role:** The 3D Sculptor.
* **Mechanism:** Using the isolated target coordinates, RFdiffusion (developed by the Baker Lab, UW) diffuses a completely new 3D protein structure around the Dengue virus target. It creates thousands of geometric "backbones" (mimetic locks) that fit the viral "keyhole".
* **Integration:** Orchestrated via `run_rfdiffusion_wrapper.py`.

### 3. Sequence Translation (ProteinMPNN)
* **The Role:** The Chemical Compiler.
* **Mechanism:** The geometric backbones from RFdiffusion lack chemical identity. ProteinMPNN reads these 3D shapes and predicts the exact sequence of amino acids required to fold into that shape stably in a biological environment (solubility, atomic forces).
* **Integration:** Handled by `parse_to_mpnn.py`.

### 4. In-Silico Validation (AlphaFold2 / ColabFold)
* **The Role:** The Uncompromising Verifier.
* **Mechanism:** The newly generated amino acid sequences are fed "blindly" into AlphaFold2. If AlphaFold computationally folds the sequence and predicts that it forcefully binds to the Dengue target motif, the molecule passes validation.

---

## 🚀 Environment Setup & Hardware

This pipeline requires a dedicated Linux environment with CUDA support to interface with PyTorch and the heavy neural network weights.

**Prerequisites (Windows Users):**
1. Ensure hardware virtualization is enabled in your BIOS.
2. Install the Windows Subsystem for Linux (WSL2):
   ```powershell
   wsl --install -d Ubuntu
   ```
3. Inside WSL2, install the NVIDIA CUDA Toolkit and PyTorch.
4. Clone this repository into the Linux subsystem to avoid I/O bottlenecks.

## 📥 External Model Weights
*Note: The actual neural network intelligence is not stored in this repository.*
Before running the generation phase, you must download the respective pretrained model weights for:
*   **RFdiffusion** (~5 GB)
*   **ProteinMPNN** (~200 MB)
*   **ColabFold/AlphaFold2** (Dependent on MSA databases)

---

## ⚖️ License
This project is licensed under the **MIT License**. 

We believe that computational approaches to severe global health threats should remain open-source. You are free to use, modify, distribute, and build upon this architecture without restriction, provided that the original copyright notice is included.

See the `LICENSE` file for full details.
