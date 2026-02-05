#!/usr/bin/env python3
"""
confidence_ambiguity_qwen_audit.py

Deterministic confidence scoring + ambiguity handling + Qwen audit for cluster cell-type prediction.

WHAT IT DOES
------------
For each cluster:
  (A) Builds candidate label scores from marker genes using a gene->cell_type knowledge base
      (e.g., PanglaoDB / CellMarker-derived JSON).
  (B) Applies ambiguity handling using:
        - mapping coverage (how many marker genes map to any cell type),
        - margin (top1 vs top2 separation),
        - entropy (how diffuse the label score distribution is).
  (C) (Optional) Audits Qwen’s predicted label using:
        - Panglao-based support score and (canonical) rank,
        - GO-slim lineage consistency score,
        - combined correctness score.

IMPORTANT FIX INCLUDED
---------------------
Qwen labels often differ from Panglao candidate labels by:
  - pluralization (e.g., "T cells" vs "T Cell")
  - parenthetical qualifiers (e.g., "Cholangiocyte (Reactive/EMT-like)")
  - minor punctuation/spacing differences

This script canonicalizes both candidate labels and Qwen labels via label_key()
so Qwen labels will be matched robustly.

INPUTS
------
Required:
  1) --markers_csv
     CSV with at least: cluster, gene, log2FC
     Recommended columns (if available): p_adj, pct_in, pct_out, cluster_size

  2) --gene_knowledge_json
     JSON mapping gene -> info, containing at least:
       gene_knowledge[GENE]["cell_type"]  (string label)
     Optional:
       gene_knowledge[GENE]["confidence_score"] (0-100)

Optional:
  3) --goslim_csv
     CSV with: cluster and a term column among:
       function | program | term_name | goslim_term

  4) --qwen_csv
     CSV with: cluster, predicted_label

OUTPUT
------
A CSV (via --out_csv) with per-cluster:
  - ambiguity-aware label and confidence
  - candidates list (JSON)
  - Qwen label + audit columns (if --qwen_csv provided)

USAGE EXAMPLE
-------------
python confidence_ambiguity_qwen_audit.py \
  --markers_csv top_20_marker_gene_per_cluster_hsnmf.csv \
  --gene_knowledge_json gene_knowledge_clean_phase1.json \
  --goslim_csv cluster_functional_programs_goslim.csv \
  --qwen_csv qwen_predictions.csv \
  --out_csv cluster_confidence_ambiguity_audit.csv
"""

from __future__ import annotations

import argparse
import json
import math
import re
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ----------------------------
# Basic utilities
# ----------------------------

def safe_float(x, default=None):
    try:
        if x is None:
            return default
        if isinstance(x, float) and math.isnan(x):
            return default
        return float(x)
    except Exception:
        return default


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))


def clamp01(x: float) -> float:
    return clamp(x, 0.0, 1.0)


def normalize_gene(g: str) -> str:
    return str(g).strip().upper()


def normalize_label(s: str) -> str:
    if s is None:
        return ""
    s = str(s).strip()
    s = re.sub(r"\s+", " ", s)
    return s


def softmax(scores: List[float], temp: float = 1.0) -> List[float]:
    if not scores:
        return []
    t = max(1e-6, float(temp))
    m = max(scores)
    exps = [math.exp((s - m) / t) for s in scores]
    z = sum(exps)
    return [e / z for e in exps]


def shannon_entropy(probs: List[float]) -> float:
    eps = 1e-12
    s = 0.0
    for p in probs:
        p = max(eps, float(p))
        s -= p * math.log(p)
    return s


# ----------------------------
# Label canonicalization (FIX)
# ----------------------------

def label_key(s: str) -> str:
    """
    Canonical key for matching Qwen labels to Panglao candidate labels.
    Makes matching robust to pluralization and parenthetical qualifiers.
    """
    s = normalize_label(s).lower()

    # drop parenthetical qualifiers (e.g. "(Reactive/EMT-like)")
    s = re.sub(r"\(.*?\)", "", s).strip()

    # remove non-alphanumeric (keep spaces)
    s = re.sub(r"[^a-z0-9\s]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()

    # naive singularization for common plurals
    if s.endswith(" cells"):
        s = s[:-1]  # "cells" -> "cell"
    elif s.endswith("s") and len(s) > 3:
        s = s[:-1]

    # optional: small alias map (extend as needed)
    alias = {
        "t cell": "t cell",
        "t cells": "t cell",
        "macrophage monocyte": "macrophage",
        "macrophage / monocyte": "macrophage",
        "macrophage": "macrophage",
        "monocyte": "monocyte",
        "hepatocyte": "hepatocyte",
        "cholangiocyte": "cholangiocyte",
        "endothelial cell": "endothelial cell",
        "fibroblast": "fibroblast",
        "stromal cell": "stromal cell",
        "stroma": "stromal cell",
        "mesenchymal progenitor": "mesenchymal progenitor",
        "mesenchymal progenitors": "mesenchymal progenitor",
    }
    return alias.get(s, s)


# ----------------------------
# Marker evidence weighting
# ----------------------------

def marker_weight(
    log2fc: Optional[float],
    pct_in: Optional[float],
    pct_out: Optional[float],
    p_adj: Optional[float],
) -> float:
    """
    Deterministic marker weight in roughly [0, ~3.25], combining:
      - effect size (log2FC)
      - specificity (pct_in - pct_out)
      - significance (-log10(p_adj))

    Works if pct_in/pct_out/p_adj are missing (uses log2FC only).
    """
    w = 0.0

    # log2FC (dominant term)
    if log2fc is not None:
        w += 1.5 * clamp01(float(log2fc) / 2.5)  # saturate around ~2.5

    # prevalence gap (specificity)
    if pct_in is not None and pct_out is not None:
        gap = max(0.0, float(pct_in) - float(pct_out))
        w += 1.0 * clamp01(gap / 0.5)  # saturate around 0.5 gap

    # adjusted p-value (significance)
    if p_adj is not None and float(p_adj) > 0.0:
        s = min(6.0, max(0.0, -math.log10(float(p_adj))))  # cap at 1e-6
        w += 0.75 * (s / 6.0)

    return w


# ----------------------------
# GO-slim parsing + lineage signal
# ----------------------------

def load_goslim_terms(goslim_csv: str) -> Dict[int, List[str]]:
    """
    Loads a CSV and extracts per-cluster GO-slim term strings.
    It looks for a term column among: program, function, term_name, goslim_term
    """
    df = pd.read_csv(goslim_csv)
    df.columns = [c.strip() for c in df.columns]

    if "cluster" not in df.columns:
        raise ValueError("goslim_csv must contain a 'cluster' column")

    term_col = None
    for c in ["program", "function", "term_name", "goslim_term"]:
        if c in df.columns:
            term_col = c
            break

    if term_col is None:
        raise ValueError(
            "goslim_csv must contain one of: program | function | term_name | goslim_term"
        )

    out: Dict[int, List[str]] = {}
    for _, r in df.iterrows():
        cid = int(r["cluster"])
        t = str(r[term_col])
        if t and t.lower() != "nan":
            out.setdefault(cid, []).append(t)
    return out


def infer_lineage_from_label(label: str) -> str:
    """
    Coarse lineage mapping used only for GO-slim consistency.
    """
    l = label_key(label)  # canonical form (lowercase)

    # immune
    if any(k in l for k in ["t cell", "b cell", "plasma", "nk", "natural killer", "macrophage", "monocyte", "myeloid", "immune"]):
        return "immune"

    # metabolic (hepatocyte-like)
    if "hepatocyte" in l:
        return "metabolic"

    # epithelial
    if any(k in l for k in ["cholangiocyte", "duct", "epithel"]):
        return "epithelial"

    # stromal
    if any(k in l for k in ["stromal", "stroma", "fibro", "mesench", "stellate", "pericyte", "endothelial"]):
        return "stromal"

    return "other"


def goslim_signal_scores(terms: Optional[List[str]]) -> Dict[str, float]:
    """
    Converts GO-slim terms into a coarse lineage signal distribution.
    Simple keyword tallying (not a program rule engine).
    """
    if not terms:
        return {"immune": 0.0, "metabolic": 0.0, "epithelial": 0.0, "stromal": 0.0, "other": 1.0}

    txt = " | ".join([str(t).lower() for t in terms])

    immune = 0.0
    metabolic = 0.0
    epithelial = 0.0
    stromal = 0.0

    for k in ["immune", "leukocyte", "t cell", "lymphocyte", "antigen", "cytokine", "chemokine", "inflammatory", "defense response", "innate", "adaptive"]:
        if k in txt:
            immune += 1.0

    for k in ["lipid", "fatty acid", "cholesterol", "gluconeogenesis", "bile", "xenobiotic", "drug metabolic", "steroid", "retinol", "alcohol metabolic", "coagulation", "hemostasis"]:
        if k in txt:
            metabolic += 1.0

    for k in ["epithelial", "cell adhesion", "cell-cell adhesion", "junction", "tight junction", "epithelium"]:
        if k in txt:
            epithelial += 1.0

    for k in ["extracellular matrix", "matrix organization", "collagen", "angiogenesis", "wound healing", "fibroblast"]:
        if k in txt:
            stromal += 1.0

    tot = immune + metabolic + epithelial + stromal
    if tot <= 0:
        return {"immune": 0.0, "metabolic": 0.0, "epithelial": 0.0, "stromal": 0.0, "other": 1.0}

    return {
        "immune": immune / tot,
        "metabolic": metabolic / tot,
        "epithelial": epithelial / tot,
        "stromal": stromal / tot,
        "other": 0.0,
    }


# ----------------------------
# Candidate scoring from gene knowledge
# ----------------------------

def compute_candidate_scores(
    markers_df: pd.DataFrame,
    gene_knowledge: Dict[str, dict],
    cluster_id: int,
    goslim_terms: Optional[List[str]] = None,
) -> Tuple[Dict[str, float], Dict[str, List[str]], float]:
    """
    Returns:
      - scores[label_key] -> float (higher is better)
      - support_genes[label_key] -> list of marker genes contributing
      - coverage -> fraction of markers mapped to any label via gene_knowledge

    NOTE:
      We store scores under canonical label_key to ensure consistency with Qwen labels.
      We also track a representative original label string separately in main().
    """
    sub = markers_df[markers_df["cluster"] == cluster_id].copy()
    if sub.empty:
        return {}, {}, 0.0

    scores: Dict[str, float] = {}
    support_genes: Dict[str, List[str]] = {}
    mapped = 0

    # Optional: tiny GO-slim uniform boost when GO evidence is strong
    go_bonus = 0.0
    if goslim_terms:
        sig = goslim_signal_scores(goslim_terms)
        go_bonus = 0.05 * max(sig["immune"], sig["metabolic"], sig["epithelial"], sig["stromal"])

    for _, row in sub.iterrows():
        gene = normalize_gene(row["gene"])
        rec = gene_knowledge.get(gene)
        if not rec:
            continue

        raw_label = normalize_label(rec.get("cell_type", ""))
        if not raw_label:
            continue

        mapped += 1

        # marker evidence weight
        w = marker_weight(
            log2fc=safe_float(row.get("log2FC")),
            pct_in=safe_float(row.get("pct_in")),
            pct_out=safe_float(row.get("pct_out")),
            p_adj=safe_float(row.get("p_adj")),
        )

        # DB confidence if provided (typical 85..100); default 100
        db_conf = safe_float(rec.get("confidence_score"), 100.0)
        db_w = clamp01((db_conf - 80.0) / 20.0)  # 80->0, 100->1

        contrib = w * (0.6 + 0.4 * db_w)

        k = label_key(raw_label)
        scores[k] = scores.get(k, 0.0) + contrib
        support_genes.setdefault(k, []).append(gene)

    coverage = mapped / max(1, len(sub))

    if go_bonus > 0 and scores:
        for k in list(scores.keys()):
            scores[k] *= (1.0 + go_bonus)

    return scores, support_genes, coverage


# ----------------------------
# Ambiguity handling + confidence
# ----------------------------

def decide_label_from_scores(
    scores: Dict[str, float],
    support_genes: Dict[str, List[str]],
    coverage: float,
    *,
    min_coverage: float,
    min_top_score: float,
    min_margin: float,
    max_entropy: float,
) -> Tuple[str, float, str, dict]:
    """
    Returns:
      predicted_label_key (or "ambiguous"/"unresolved"),
      confidence in [0,1],
      status in {HIGH, MEDIUM, AMBIGUOUS, UNRESOLVED},
      diagnostics dict
    """
    if not scores or coverage < min_coverage:
        return "unresolved", 0.0, "UNRESOLVED", {
            "reason": "insufficient_mapping_coverage",
            "coverage": float(coverage),
        }

    ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
    top1_key, top1_score = ranked[0]
    top2_key, top2_score = ranked[1] if len(ranked) > 1 else ("", 0.0)

    probs = softmax([s for _, s in ranked], temp=1.0)
    ent = shannon_entropy(probs)

    denom = max(1e-6, float(top1_score))
    margin = (float(top1_score) - float(top2_score)) / denom

    # confidence scalar
    strength = clamp01(float(top1_score) / 4.0)  # saturate around 4
    conf = 0.45 * strength + 0.25 * clamp01(margin / 0.5) + 0.20 * clamp01((coverage - min_coverage) / (1.0 - min_coverage))
    conf *= (1.0 - 0.25 * clamp01((ent - 0.5) / 1.5))
    conf = clamp01(conf)

    sup1 = len(set(support_genes.get(top1_key, [])))

    if float(top1_score) < min_top_score:
        return "ambiguous", conf, "AMBIGUOUS", {
            "reason": "low_top_score",
            "top1_key": top1_key,
            "top1_score": float(top1_score),
            "top2_key": top2_key,
            "top2_score": float(top2_score),
            "coverage": float(coverage),
            "entropy": float(ent),
            "margin": float(margin),
            "supporting_markers_top1": sup1,
        }

    if margin < min_margin:
        return "ambiguous", conf, "AMBIGUOUS", {
            "reason": "small_margin",
            "top1_key": top1_key,
            "top1_score": float(top1_score),
            "top2_key": top2_key,
            "top2_score": float(top2_score),
            "coverage": float(coverage),
            "entropy": float(ent),
            "margin": float(margin),
            "supporting_markers_top1": sup1,
        }

    if ent > max_entropy:
        return "ambiguous", conf, "AMBIGUOUS", {
            "reason": "high_entropy",
            "top1_key": top1_key,
            "top1_score": float(top1_score),
            "top2_key": top2_key,
            "top2_score": float(top2_score),
            "coverage": float(coverage),
            "entropy": float(ent),
            "margin": float(margin),
            "supporting_markers_top1": sup1,
        }

    status = "HIGH" if conf >= 0.70 else "MEDIUM"
    return top1_key, conf, status, {
        "reason": "accepted",
        "top1_key": top1_key,
        "top1_score": float(top1_score),
        "top2_key": top2_key,
        "top2_score": float(top2_score),
        "coverage": float(coverage),
        "entropy": float(ent),
        "margin": float(margin),
        "supporting_markers_top1": sup1,
    }


# ----------------------------
# Qwen audit (correctness)
# ----------------------------

def evaluate_qwen_label(
    qwen_label: str,
    candidate_scores: Dict[str, float],      # keyed by label_key
    goslim_terms: Optional[List[str]],
) -> Dict[str, object]:
    """
    Uses canonical label keys for robust matching.
    """
    qlab = normalize_label(qwen_label)
    qk = label_key(qlab) if qlab else ""

    if not qk:
        return {
            "panglao_support": 0.0,
            "panglao_rank": -1,
            "panglao_margin": 0.0,
            "goslim_lineage": "other",
            "goslim_consistency": 0.0,
            "combined_correctness": 0.0,
            "flag": "missing_qwen_label",
        }

    ranked = sorted(candidate_scores.items(), key=lambda x: x[1], reverse=True)
    if not ranked:
        pang_support = 0.0
        pang_rank = -1
        pang_margin = 0.0
        best_key = ""
        best_score = 0.0
    else:
        best_key, best_score = ranked[0]

        q_score = float(candidate_scores.get(qk, 0.0))
        pang_support = clamp01(q_score / 4.0)  # saturate around ~4

        # rank among canonical keys
        keys_sorted = [k for k, _ in ranked]
        pang_rank = (1 + keys_sorted.index(qk)) if qk in candidate_scores else -1

        pang_margin = (q_score - float(best_score)) / float(best_score) if float(best_score) > 1e-9 else 0.0

    lineage = infer_lineage_from_label(qlab)
    sig = goslim_signal_scores(goslim_terms)
    gos_cons = float(sig.get(lineage, 0.0))

    combined = 0.75 * float(pang_support) + 0.25 * float(gos_cons)

    flag = "ok"
    if pang_rank == -1:
        flag = "not_supported_by_panglao_candidates"
    elif pang_support < 0.35:
        flag = "weak_panglao_support"
    elif lineage != "other" and gos_cons < 0.20:
        flag = "goslim_inconsistent"

    return {
        "panglao_support": round(float(pang_support), 4),
        "panglao_rank": int(pang_rank),
        "panglao_margin": round(float(pang_margin), 4),
        "goslim_lineage": lineage,
        "goslim_consistency": round(float(gos_cons), 4),
        "combined_correctness": round(float(combined), 4),
        "flag": flag,
        "best_candidate_key": best_key,
        "best_candidate_score": round(float(best_score), 4),
    }


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--markers_csv", required=True, help="Top-K markers per cluster CSV")
    ap.add_argument("--gene_knowledge_json", required=True, help="Gene->cell_type knowledge JSON (Panglao/CellMarker)")
    ap.add_argument("--goslim_csv", default=None, help="Optional GO-slim per-cluster CSV")
    ap.add_argument("--qwen_csv", default=None, help="Optional Qwen predictions CSV: cluster,predicted_label")
    ap.add_argument("--out_csv", required=True, help="Output CSV path")

    # ambiguity thresholds
    ap.add_argument("--min_coverage", type=float, default=0.25)
    ap.add_argument("--min_top_score", type=float, default=1.25)
    ap.add_argument("--min_margin", type=float, default=0.20)
    ap.add_argument("--max_entropy", type=float, default=1.35)

    args = ap.parse_args()

    # markers
    markers = pd.read_csv(args.markers_csv)
    markers.columns = [c.strip() for c in markers.columns]
    required_cols = {"cluster", "gene", "log2FC"}
    missing = required_cols - set(markers.columns)
    if missing:
        raise ValueError(f"markers_csv missing required columns: {sorted(list(missing))}")

    markers["cluster"] = markers["cluster"].astype(int)
    markers["gene"] = markers["gene"].astype(str).apply(normalize_gene)

    # knowledge
    with open(args.gene_knowledge_json, "r") as f:
        gene_knowledge_raw = json.load(f)
    gene_knowledge = {normalize_gene(k): v for k, v in gene_knowledge_raw.items()}

    # GO-slim terms (optional)
    goslim_map: Dict[int, List[str]] = {}
    if args.goslim_csv:
        goslim_map = load_goslim_terms(args.goslim_csv)

    # Qwen predictions (optional)
    qwen_df = None
    if args.qwen_csv:
        qwen_df = pd.read_csv(args.qwen_csv)
        qwen_df.columns = [c.strip() for c in qwen_df.columns]
        if "cluster" not in qwen_df.columns or "predicted_label" not in qwen_df.columns:
            raise ValueError("qwen_csv must contain columns: cluster,predicted_label")
        qwen_df["cluster"] = qwen_df["cluster"].astype(int)
        qwen_df["predicted_label"] = qwen_df["predicted_label"].astype(str)

    clusters = sorted(markers["cluster"].unique().tolist())
    rows: List[Dict[str, object]] = []

    for cid in clusters:
        terms = goslim_map.get(cid)

        # Candidate scores keyed by canonical label_key
        scores, support, cov = compute_candidate_scores(
            markers_df=markers,
            gene_knowledge=gene_knowledge,
            cluster_id=cid,
            goslim_terms=terms,
        )

        # Representative labels for display: pick the most common raw label for each key
        # (best-effort; if you want perfect mapping, store raw label in gene_knowledge per key)
        key_to_display = {k: k for k in scores.keys()}  # default: show canonical key

        pred_key, conf, status, diag = decide_label_from_scores(
            scores=scores,
            support_genes=support,
            coverage=cov,
            min_coverage=args.min_coverage,
            min_top_score=args.min_top_score,
            min_margin=args.min_margin,
            max_entropy=args.max_entropy,
        )

        ranked = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        top1 = ranked[0] if ranked else ("", 0.0)
        top2 = ranked[1] if len(ranked) > 1 else ("", 0.0)

        # Candidate list (top 10) with support genes
        candidate_list = [
            {
                "label_key": k,
                "score": round(float(sc), 6),
                "support_genes": sorted(list(set(support.get(k, [])))),
            }
            for k, sc in ranked[:10]
        ]

        # Cluster size (if present)
        cluster_size = None
        if "cluster_size" in markers.columns:
            try:
                cluster_size = int(markers[markers["cluster"] == cid]["cluster_size"].iloc[0])
            except Exception:
                cluster_size = None

        # Qwen label (optional)
        qwen_label = ""
        if qwen_df is not None:
            hit = qwen_df[qwen_df["cluster"] == cid]
            if len(hit) > 0:
                qwen_label = str(hit.iloc[0]["predicted_label"])

        audit = {}
        if qwen_df is not None:
            audit = evaluate_qwen_label(qwen_label, scores, terms)

        # Convert internal keys to human-friendly final label strings
        if pred_key == "ambiguous":
            pred_label = "Ambiguous"
        elif pred_key == "unresolved":
            pred_label = "Unresolved"
        else:
            pred_label = pred_key  # canonical key is already readable ("t cell", "hepatocyte", etc.)

        row = {
            "cluster": cid,
            "cluster_size": cluster_size,
            "predicted_label_after_ambiguity": pred_label,
            "confidence": round(float(conf), 4),
            "status": status,
            "top1_label_key": top1[0],
            "top1_score": round(float(top1[1]), 6),
            "top2_label_key": top2[0],
            "top2_score": round(float(top2[1]), 6),
            "marker_coverage": round(float(diag.get("coverage", cov)), 4),
            "margin": round(float(diag.get("margin", 0.0)), 4),
            "entropy": round(float(diag.get("entropy", 0.0)), 4),
            "supporting_markers_top1": int(diag.get("supporting_markers_top1", 0)),
            "diagnostics_json": json.dumps(diag),
            "candidates_json": json.dumps(candidate_list),
        }

        if qwen_df is not None:
            row.update({
                "qwen_label": qwen_label,
                "qwen_label_key": label_key(qwen_label) if qwen_label else "",
                **audit,
            })

        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(args.out_csv, index=False)
    print(f"Wrote: {args.out_csv}")


if __name__ == "__main__":
    main()
