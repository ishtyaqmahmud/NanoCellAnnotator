# NanoCellAnnotator

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
