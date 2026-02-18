# NanoCellAnnotator: Formalizing Expert Cell Type Annotation with Large Language Models

**NanoCellAnnotator** is a biologically constrained, confidence-aware framework for automated cell-type annotation in spatial transcriptomics. Unlike unconstrained LLM approaches that risk biologically unsupported predictions, NanoCellAnnotator decouples spatial discovery from semantic inference to ensure deterministic, interpretable, and reproducible results.


## 🔬 Core Innovations

- **Decoupled Architecture**: NanoCellAnnotator explicitly separates spatial structure discovery (via hSNMF) from semantic labeling to ensure annotations are grounded in spatially coherent cellular domains.

- **Ontology-Grounded Reasoning**: The framework projects cluster-specific marker genes into functional programs using Gene Ontology (GO) enrichment and GO-slim projection, providing structured biological context without premature cell-type assignment.

- **Constrained Semantic Inference**: The system eliminates LLM hallucinations by restricting the admissible label space to curated entries from PanglaoDB and CellMarker.

- **Independent Uncertainty Quantification**: A deterministic confidence scoring mechanism evaluates marker-gene support and lineage separation independently of the language model to explicitly flag ambiguous clusters.

- **Edge-Ready Deployment**: The pipeline utilizes Qwen2.5-1.5B-Instruct, a lightweight, locally executable model that ensures reproducibility and data privacy on commodity-grade hardware.


## ⚙️ Methodology Summary

NanoCellAnnotator follows a deterministic, multi-stage pipeline:

1. **Upstream Spatial Structure Discovery**
   Spatially coherent cellular domains are identified using Hybrid Spatially Regularized Non-negative Matrix Factorization (hSNMF). This stage is model-agnostic and operates independently of any cell-type labels to define geometrically contiguous domains.

2. **Deterministic Biological Evidence Construction**
   Cluster-specific marker genes are extracted via differential expression analysis. These markers are mapped to ontology-derived functional programs using GO enrichment and GO-slim projection to provide stable summaries of biological activity.

3. **Constrained Semantic Inference**
   A lightweight, locally executed LLM (Qwen2.5-1.5B) selects a single cell-type label per cluster. The inference is strictly bounded by an admissible label space derived from PanglaoDB and CellMarker to prevent biologically unsupported predictions.

4. **Confidence-Aware Annotation**
   Annotation confidence is assessed independently of the LLM by quantifying marker-gene support strength and lineage separation. This allows the framework to explicitly flag ambiguous or heterogeneous clusters as "Unresolved" rather than forcing a label.

## 🧬 Framework Overview

<img width="909" height="452" alt="Methodology" src="https://github.com/user-attachments/assets/9cef378d-424b-4b8f-81b7-5045a7574ff8" />

**Figure:** Schematic of the NanoCellAnnotator framework. The pipeline decouples spatial clustering from semantic inference to ensure reproducibility. (1) Spatial Structure Discovery: Raw data is processed via hybrid spatially regularized NMF (hSNMF) to identify coherent tissue domains. (2) Biological Evidence Construction: Marker genes are mapped to functional GO terms and constrained by reference databases (PanglaoDB, CellMarker). (3) Constrained Semantic Inference: A lightweight LLM (NanoLLM) assigns labels using only structured summaries—without access to raw expression—to ensure deterministic outputs. (4) Confidence Assessment: Final annotations are validated against marker support to explicitly flag ambiguous clusters.

## 📋 Data Availability

The clinical Xenium spatial transcriptomics datasets used in this study were obtained from MD Anderson Cancer Center and are subject to institutional and patient privacy restrictions. Due to these constraints, raw clinical data cannot be publicly released.

To ensure transparency and reproducibility, this repository provides:

- The complete NanoCellAnnotator pipeline implementation
- All curated public reference resources (PanglaoDB, CellMarker, GO-slim)

Researchers with access to appropriate Xenium datasets may apply NanoCellAnnotator directly to their own data following the provided instructions.

## 🏛️ Repository Navigation Guide

The repository is organized by research stages rather than as a traditional software package.

### 1. Public Biological Reference Resources (`Data/`)
Contains curated, public biological knowledge used to constrain inference:
- PanglaoDB marker gene database
- CellMarker reference annotations
- Gene Ontology (GO) and GO-slim resources

These resources are used exclusively to restrict the admissible cell-type label space and do not contain clinical data.

---

### 2. Evidence Construction and Knowledge Processing (`Data_processing/`)
Implements deterministic preprocessing steps including:
- Marker gene processing
- Ontology-based functional abstraction
- Construction of structured biological summaries

This stage corresponds to Sections 2.2–2.3 of the paper.

---

### 3. Constrained Inference and Annotation (`model_run/`, `Phase2/`)
Implements the core NanoLLM-guided annotation logic:
- Candidate label restriction
- Deterministic label selection
- Integration of marker and functional evidence

These components correspond to Section 2.4 of the paper.

---

### 4. Confidence Scoring and Ambiguity Handling (`confidence/`)
Implements independent, deterministic confidence assessment:
- Marker support aggregation
- Lineage separation metrics
- Ambiguity flagging

This stage is independent of LLM inference and corresponds to Section 2.5 of the paper.

---

### 5. Ablation Studies (`Abluation_study/`)
Contains notebooks and outputs used to generate ablation results reported in the paper (Table 2).

These experiments isolate the effects of evidence representation and constraint mechanisms.

## 🚀 Installation

This section provides step-by-step instructions to set up **NanoCellAnnotator** in a clean and reproducible Python environment. The implementation is designed to be platform-agnostic and has been tested on macOS and Linux systems.

---

### 1. Clone the Repository

Begin by downloading the NanoCellAnnotator codebase from GitHub:

```bash
git clone [https://github.com/ishtyaqmahmud/NanoCellAnnotator.git](https://github.com/ishtyaqmahmud/NanoCellAnnotator.git)
cd NanoCellAnnotator
```
### 2. Create a Virtual Environment 
To avoid dependency conflicts, we strongly recommend using a dedicated Python 3.10 environment.

```bash
conda create -n nano_llm_env python=3.10 -y
conda activate nano_llm_env
```

### 3. Install Dependencies
All required high-level Python dependencies are listed in the requirements.txt file located at the repository root. Install them via pip:
```bash
pip install --upgrade pip
pip install -r requirements.txt
```


## Reference

If you find NanoCellAnnotator useful in your research, please consider citing our paper:

Mahmud, M.I., Anzum, H., Kochat, V., Satpati, S., Rai, K., and Banerjee, T. *NanoCellAnnotator: Formalizing Expert Cell Type Annotation with Large Language Models*. ISMB 2026 (submitted).

```bibtex
@misc{mahmud2026nanocellannotator,
  title = {NanoCellAnnotator: Formalizing Expert Cell Type Annotation with Large Language Models},
  author = {Mahmud, M. I. and Anzum, H. and Kochat, V. and Satpati, S. and Rai, K. and Banerjee, T.},
  year = {2026},
  note = {Submitted to ISMB 2026}
}
```
