# NanoCellAnnotator: Formalizing Expert Cell Type Annotation with Large Language Models

## Overview

NanoCellAnnotator is a biologically constrained and confidence-aware framework for automated cell-type annotation in spatial transcriptomics data. The method is designed to address key challenges in spatial analysis, including data sparsity, spatial heterogeneity, and the absence of comprehensive labeled reference atlases.

Unlike unconstrained LLM-based approaches, NanoCellAnnotator explicitly decouples spatial structure discovery, deterministic biological evidence construction, and language-model–based semantic inference into separate stages. Spatially coherent cellular domains are identified upstream, while cluster-level marker genes are abstracted into ontology-derived functional programs using Gene Ontology enrichment and GO-slim projection.

A lightweight large language model (NanoLLM) is employed solely as a constrained semantic integrator, operating on structured biological summaries rather than raw gene expression data or spatial coordinates. Curated reference databases (PanglaoDB and CellMarker) are used exclusively to restrict the admissible space of cell-type labels, preventing biologically unsupported predictions.

Importantly, annotation confidence is assessed independently of language-model inference based on marker-gene support strength and lineage separation. This allows ambiguous or heterogeneous spatial clusters to be explicitly flagged rather than forcibly labeled, preserving biologically meaningful uncertainty.

NanoCellAnnotator operates in a fully unsupervised setting and does not rely on labeled training data or reference atlases. The framework emphasizes reproducibility, interpretability, and evidence-consistent reasoning, making it suitable for both academic research and translational applications in spatial transcriptomics.


## Methodology Summary

NanoCellAnnotator follows a deterministic, multi-stage pipeline:

1. **Upstream Spatial Structure Discovery**
   Spatially coherent cellular domains are identified using spatially regularized clustering, independent of any cell-type labels.

2. **Deterministic Biological Evidence Construction**
   Cluster-specific marker genes are extracted and mapped to ontology-derived functional programs using Gene Ontology enrichment and GO-slim projection. Curated databases are used only to restrict admissible cell-type labels.

3. **Constrained Semantic Inference**
   A lightweight LLM (NanoLLM) selects a single cell-type label per cluster from the restricted label space based solely on structured biological summaries.

4. **Confidence-Aware Annotation**
   Marker-gene support strength and lineage separation are evaluated independently to determine annotation confidence and flag ambiguous clusters.

## Framework Overview

<img width="909" height="452" alt="Methodology" src="https://github.com/user-attachments/assets/9cef378d-424b-4b8f-81b7-5045a7574ff8" />

**Figure:** Schematic of the NanoCellAnnotator framework. The pipeline decouples spatial clustering from semantic inference to ensure reproducibility. (1) Spatial Structure Discovery: Raw data is processed via hybrid spatially regularized NMF (hSNMF) to identify coherent tissue domains. (2) Biological Evidence Construction: Marker genes are mapped to functional GO terms and constrained by reference databases (PanglaoDB, CellMarker). (3) Constrained Semantic Inference: A lightweight LLM (NanoLLM) assigns labels using only structured summaries—without access to raw expression—to ensure deterministic outputs. (4) Confidence Assessment: Final annotations are validated against marker support to explicitly flag ambiguous clusters.

## Data Availability

The clinical Xenium spatial transcriptomics datasets used in this study were obtained from MD Anderson Cancer Center and are subject to institutional and patient privacy restrictions. Due to these constraints, raw clinical data cannot be publicly released.

To ensure transparency and reproducibility, this repository provides:

- The complete NanoCellAnnotator pipeline implementation
- All curated public reference resources (PanglaoDB, CellMarker, GO-slim)

Researchers with access to appropriate Xenium datasets may apply NanoCellAnnotator directly to their own data following the provided instructions.

## Repository Navigation Guide

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

If you use NanoCellAnnotator, please cite:

Mahmud, M.I., Anzum, H., Kochat, V., Satpati, S., Rai, K., and Banerjee, T. *NanoCellAnnotator: Formalizing Expert Cell Type Annotation with Large Language Models*. ISMB 2026 (submitted).
