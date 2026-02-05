#!/usr/bin/env python3
"""
Updated evaluation script.

Changes made:
- Filters out generic 'normal cell'-like candidates before ranking (as before).
- Computes Top-k agreement (cell-weighted) on:
    * All clusters
    * High + Moderate confidence clusters (primary reported metric)
    * High-only clusters (sensitivity)
  If no explicit confidence tier column is present, the script will (1) try to
  derive tiers from a numeric "confidence" column using conservative defaults,
  and (2) otherwise fall back to evaluating on all clusters while warning the user.
- Prints commitment rates (fraction of cells retained) for each mask.
- Saves per-cluster audit CSV with top-k membership flags.
- Robust parsing of candidates_json and tolerant to missing fields.
"""

import pandas as pd
import json
import ast
import sys
from typing import List

# =========================
# Configuration
# =========================
INPUT_CSV = "out/confidence_a4.csv"            # Change as needed
OUTPUT_AUDIT_CSV = "out/topk_audit_A4_filtered_updated.csv"

# Generic / catch-all labels to exclude from candidate ranking
GENERIC_LABELS = {
    "normal cell",
    "cell",
    "other",
    "unknown",
    "uncharacterized",
    "unspecified",
}

TOPK_LIST = [1, 2, 3, 5]
N_DEBUG = 10  # number of rows to print for quick sanity checks

# Confidence thresholds (used if `confidence` numeric column exists and no tier column)
CONF_HIGH = 0.70
CONF_MODERATE = 0.50


# =========================
# Helpers
# =========================
def parse_candidates(s):
    """Parse candidates_json into a Python list of dicts."""
    if pd.isna(s):
        return []
    try:
        return json.loads(s)
    except Exception:
        try:
            return ast.literal_eval(s)
        except Exception:
            # last resort: try simple string cleaning
            return []


def filter_generic_candidates(candidates: List[dict]) -> List[dict]:
    """
    Remove generic / catch-all labels from candidate list.
    Keeps ordering of remaining candidates.
    """
    filtered = []
    for d in candidates:
        if not isinstance(d, dict):
            continue
        label = str(d.get("label_key", "")).lower().strip()
        if label and label not in GENERIC_LABELS:
            filtered.append(d)
    return filtered


def get_topk_keys(candidates: List[dict], k: int) -> List[str]:
    """Extract label_key for the top-k entries of the candidate list."""
    keys = []
    for d in candidates[:k]:
        if isinstance(d, dict) and "label_key" in d:
            keys.append(d["label_key"])
    return keys


def compute_cell_weighted_pct(df: pd.DataFrame, mask: pd.Series, flag_col: str) -> float:
    total_cells = df.loc[mask, "cluster_size"].sum()
    if total_cells == 0:
        return 0.0
    weighted_hits = (df.loc[mask, flag_col] * df.loc[mask, "cluster_size"]).sum()
    return float(weighted_hits) / float(total_cells) * 100.0


def detect_and_build_tiers(df: pd.DataFrame) -> pd.Series:
    """
    Return a Series of confidence tiers: 'High'/'Moderate'/'Low'.
    Priority:
      1) If 'confidence_tier' exists, use it (normalize strings).
      2) Else if numeric 'confidence' exists, use thresholds CONF_HIGH and CONF_MODERATE.
      3) Else if 'status' or 'flag' contains strings like 'High', try to use them.
      4) Otherwise return a Series of 'Unknown' (treated as all included).
    """
    col_candidates = ["confidence_tier", "confidence_label", "tier"]
    for c in col_candidates:
        if c in df.columns:
            ser = df[c].astype(str).str.strip().str.title()
            # normalize a few common variants
            ser = ser.replace({"Hi": "High", "Md": "Moderate", "Low ": "Low"})
            return ser

    if "confidence" in df.columns:
        # numeric mapping
        try:
            conf = pd.to_numeric(df["confidence"], errors="coerce")
            tiers = pd.Series(["Low"] * len(df), index=df.index)
            tiers.loc[conf >= CONF_HIGH] = "High"
            tiers.loc[(conf >= CONF_MODERATE) & (conf < CONF_HIGH)] = "Moderate"
            tiers.loc[conf.isna()] = "Unknown"
            return tiers
        except Exception:
            pass

    # fallback: check 'status'/'flag' textual columns
    for c in ["status", "flag"]:
        if c in df.columns:
            ser = df[c].astype(str).str.strip().str.title()
            # map obvious variants
            ser = ser.replace({"Ok": "High", "Good": "High"})
            # anything not mapping we keep as-is
            return ser

    # last resort: unknown for all
    return pd.Series(["Unknown"] * len(df), index=df.index)


# =========================
# Main
# =========================
def main(input_csv=INPUT_CSV, output_audit_csv=OUTPUT_AUDIT_CSV):
    df = pd.read_csv(input_csv)

    # Basic validation
    required_cols = {"cluster", "cluster_size", "qwen_label_key", "candidates_json"}
    missing = required_cols - set(df.columns)
    if missing:
        print(f"ERROR: input CSV is missing required columns: {missing}", file=sys.stderr)
        print("Required columns (at minimum): 'cluster','cluster_size','qwen_label_key','candidates_json'", file=sys.stderr)
        return

    # Parse raw candidates and filter out generic labels from ranking
    df["candidates_raw"] = df["candidates_json"].apply(parse_candidates)
    df["candidates"] = df["candidates_raw"].apply(filter_generic_candidates)

    # Recompute top1_label_key from filtered candidates (optional)
    df["top1_label_key_filtered"] = df["candidates"].apply(
        lambda c: c[0]["label_key"] if (len(c) > 0 and isinstance(c[0], dict) and "label_key" in c[0]) else None
    )

    # Debug prints (first N rows)
    print("\n=== Candidate parsing sanity check (first N rows) ===")
    for i, row in df.head(N_DEBUG).iterrows():
        raw_top5 = get_topk_keys(row["candidates_raw"], 5)
        filt_top5 = get_topk_keys(row["candidates"], 5)
        print(
            f"row={i:>3} | cluster={row['cluster']:>3} | size={row['cluster_size']:>6} | "
            f"qwen={row['qwen_label_key']:<30} | raw_top5={raw_top5} | filt_top5={filt_top5}"
        )

    # Show clusters where generic labels were present originally
    print("\n=== Clusters where generic labels were present in raw candidates (showing top5 before/after) ===")
    shown = 0
    for _, row in df.iterrows():
        raw_labels = [d.get("label_key") for d in row["candidates_raw"] if isinstance(d, dict)]
        raw_labels_lower = [str(x).lower() for x in raw_labels]
        if any(x in raw_labels_lower for x in GENERIC_LABELS):
            print(
                f"Cluster {row['cluster']:>3} | raw_top5={raw_labels[:5]} | filt_top5={get_topk_keys(row['candidates'], 5)}"
            )
            shown += 1
            if shown >= 25:  # avoid flooding stdout
                break

    # Add per-row top-k membership flags (after filtering)
    for k in TOPK_LIST:
        df[f"in_top{k}_filtered"] = df.apply(
            lambda row: bool(row["qwen_label_key"] in get_topk_keys(row["candidates"], k)),
            axis=1,
        )

    # Determine confidence tiers (High/Moderate/Low/Unknown)
    df["confidence_tier"] = detect_and_build_tiers(df)

    # Masks
    mask_all = pd.Series(True, index=df.index)
    mask_himod = df["confidence_tier"].isin(["High", "Moderate"])
    mask_high = df["confidence_tier"].isin(["High"])

    # If tiers are Unknown for all, fallback to include all and warn
    if (df["confidence_tier"] == "Unknown").all():
        print("\nWARNING: No confidence tier information detected. Evaluating on all clusters.")
        mask_himod = mask_all
        mask_high = mask_all

    # Print Top-k metrics for All, High+Moderate (primary), High-only (sensitivity)
    print("\n=== Cell-weighted Top-k agreement (filtered candidates) ===")
    for k in TOPK_LIST:
        pct_all = compute_cell_weighted_pct(df, mask_all, f"in_top{k}_filtered")
        pct_himod = compute_cell_weighted_pct(df, mask_himod, f"in_top{k}_filtered")
        pct_high = compute_cell_weighted_pct(df, mask_high, f"in_top{k}_filtered")

        # Commitment rates
        commit_himod = float(df.loc[mask_himod, "cluster_size"].sum()) / float(df["cluster_size"].sum()) * 100.0
        commit_high = float(df.loc[mask_high, "cluster_size"].sum()) / float(df["cluster_size"].sum()) * 100.0

        print(
            f"Top-{k}: All={pct_all:6.2f}% | High+Mod={pct_himod:6.2f}% (commit={commit_himod:.1f}%) | HighOnly={pct_high:6.2f}% (commit={commit_high:.1f}%)"
        )

    # Save a cluster-level audit table
    audit_cols = [
        "cluster",
        "cluster_size",
        "qwen_label_key",
        "top1_label_key",            # original (if present)
        "top1_label_key_filtered",   # recomputed after filtering
        "confidence_tier",
    ]

    # include top-k flags in audit
    for k in TOPK_LIST:
        audit_cols.append(f"in_top{k}_filtered")

    audit = df[audit_cols].copy()
    audit = audit.sort_values("cluster")
    audit.to_csv(output_audit_csv, index=False)
    print(f"\nSaved cluster-level top-k audit table to: {output_audit_csv}")

    # Optional: print largest misses (not in top-5 after filtering)
    print("\n=== Largest misses after filtering (Qwen not in top-5), by cluster_size ===")
    misses = audit.loc[~audit["in_top5_filtered"]].sort_values("cluster_size", ascending=False)
    if not misses.empty:
        print(misses.head(20).to_string(index=False))
    else:
        print("None (all Qwen labels are in top-5 after filtering).")


if __name__ == "__main__":
    main()

