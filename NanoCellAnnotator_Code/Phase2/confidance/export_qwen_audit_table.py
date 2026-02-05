import pandas as pd

# Load A1 results
df = pd.read_csv("out/confidence_a1.csv")

# -------- Normalization --------
def normalize_label(s):
    if pd.isna(s):
        return ""
    s = s.lower()
    for ch in ["(", ")", "/", "-", ","]:
        s = s.replace(ch, " ")
    return " ".join(s.split())

# -------- Lineage assignment --------
def assign_lineage(label):
    if label is None or label == "":
        return "other"
    label = label.lower()

    if any(x in label for x in ["t cell", "t lymphocyte", "macrophage", "monocyte", "dendritic"]):
        return "immune"
    if "hepatocyte" in label or "liver cell" in label:
        return "hepatocyte"
    if any(x in label for x in ["cholangiocyte", "ductal", "basal"]):
        return "cholangiocyte"
    if any(x in label for x in ["epithelial", "enterocyte", "airway epithelial"]):
        return "epithelial"
    if any(x in label for x in ["fibroblast", "stroma", "stromal", "mesenchymal", "adipocyte"]):
        return "stromal"
    return "other"

# Apply normalization and lineage mapping
df["qwen_norm"] = df["qwen_label_key"].apply(normalize_label)
df["top1_norm"] = df["top1_label_key"].apply(normalize_label)

df["qwen_lineage"] = df["qwen_norm"].apply(assign_lineage)
df["top1_lineage"] = df["top1_norm"].apply(assign_lineage)

# -------- Cluster-level audit table --------
audit_table = (
    df[
        [
            "cluster",
            "cluster_size",
            "qwen_label_key",
            "qwen_lineage",
            "top1_label_key",
            "top1_lineage",
        ]
    ]
    .sort_values("cluster")
)

# Save for inspection / supplement
audit_table.to_csv("qwen_cluster_audit_A1.csv", index=False)

print(audit_table)

