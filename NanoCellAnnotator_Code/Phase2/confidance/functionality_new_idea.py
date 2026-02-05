from __future__ import annotations

"""Deterministic functional-program inference for marker genes.

This script performs over-representation analysis (ORA) on marker genes against
biological gene-set libraries and emits functional programs with full
provenance.

If you provide a GO DAG (go-basic.obo) and a GO-slim subset (e.g.
goslim_generic.obo), the script projects significant GO terms to GO-slim and
uses GO-slim term names as tissue-agnostic program labels (no hardcoded rules).

If GO-slim assets are not provided, the script outputs clustered enriched terms
as labels (still with stable IDs when available).

Outputs
- cluster
- function (program label)
- supporting_genes (semicolon separated)
- supporting_terms_json (JSON list of evidence terms)

Notes
- For best GO behavior, use the GOA backend (GAF + OBO) rather than MSigDB C5.
- For MSigDB GMTs, stable IDs are preserved when present in the term/description.
"""

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


# ----------------------------
# Data structures
# ----------------------------


@dataclass(frozen=True)
class TermSet:
    term_id: str
    term_name: str
    source: str
    genes: Set[str]


@dataclass(frozen=True)
class EnrichHit:
    term_id: str
    term_name: str
    source: str
    set_size: int
    overlap: int
    overlap_genes: Set[str]
    pval: float
    fdr: float


# ----------------------------
# Helpers
# ----------------------------


_RE_GO_ID = re.compile(r"\bGO:\d{7}\b")
_RE_REACTOME_ID = re.compile(r"\bR-HSA-\d+\b", re.IGNORECASE)


def norm_gene(g: str) -> str:
    return g.strip().upper()


def infer_source_from_term(raw_term: str) -> str:
    t = raw_term.upper()
    if t.startswith("GO_BIOLOGICAL_PROCESS_") or t.startswith("GOBP_") or t.startswith("GOBP "):
        return "GO:BP"
    if t.startswith("GO_MOLECULAR_FUNCTION_") or t.startswith("GOMF_") or t.startswith("GOMF "):
        return "GO:MF"
    if t.startswith("GO_CELLULAR_COMPONENT_") or t.startswith("GOCC_") or t.startswith("GOCC "):
        return "GO:CC"
    if t.startswith("REACTOME_"):
        return "Reactome"
    if t.startswith("KEGG_"):
        return "KEGG"
    if t.startswith("HALLMARK_"):
        return "MSigDB:H"
    return "Unknown"


def prettify_term_name(raw: str) -> str:
    t = raw
    for prefix in [
        "HALLMARK_",
        "REACTOME_",
        "KEGG_",
        "GO_BIOLOGICAL_PROCESS_",
        "GO_MOLECULAR_FUNCTION_",
        "GO_CELLULAR_COMPONENT_",
        "GOBP_",
        "GOMF_",
        "GOCC_",
    ]:
        if t.upper().startswith(prefix):
            t = t[len(prefix) :]
            break
    t = t.replace("_", " ").strip()
    words = t.lower().split()
    small = {"and", "or", "of", "to", "in", "by", "with", "for", "via"}
    out = [w if w in small else w.capitalize() for w in words]
    return " ".join(out)


def stable_fallback_id(raw_term: str, source: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", raw_term.strip().upper()).strip("_")
    return f"{source}:{slug}" if source != "Unknown" else f"TERM:{slug}"


def parse_ids_from_gmt_header(term: str, desc: str) -> Optional[str]:
    m_go = _RE_GO_ID.search(term) or _RE_GO_ID.search(desc)
    if m_go:
        return m_go.group(0)
    m_re = _RE_REACTOME_ID.search(term) or _RE_REACTOME_ID.search(desc)
    if m_re:
        return m_re.group(0).upper()
    return None


# ----------------------------
# I/O: GMT parsing
# ----------------------------


def load_gmt(gmt_path: str | Path, source_override: Optional[str] = None) -> Dict[str, TermSet]:
    """Parse a GMT file and return term_id -> TermSet."""
    gmt_path = Path(gmt_path)
    termsets: Dict[str, TermSet] = {}

    with gmt_path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            raw_term = parts[0]
            desc = parts[1]
            genes = {norm_gene(g) for g in parts[2:] if g.strip()}
            if not genes:
                continue

            source = source_override or infer_source_from_term(raw_term)
            extracted_id = parse_ids_from_gmt_header(raw_term, desc)
            term_id = extracted_id or stable_fallback_id(raw_term, source)
            term_name = prettify_term_name(raw_term)

            termsets[term_id] = TermSet(term_id=term_id, term_name=term_name, source=source, genes=genes)

    return termsets


# ----------------------------
# I/O: GO OBO parsing
# ----------------------------


@dataclass
class GOTerm:
    go_id: str
    name: str
    parents: Set[str]
    depth: int = 0


def load_go_dag(go_obo_path: str | Path) -> Dict[str, GOTerm]:
    """Parse go-basic.obo (or equivalent) and return a minimal GO DAG.

    We parse only:
      - id
      - name
      - is_a parents

    This is sufficient for GO-slim projection and ancestor traversal.
    """
    go_obo_path = Path(go_obo_path)

    dag: Dict[str, GOTerm] = {}
    current_id: Optional[str] = None
    current_name: Optional[str] = None
    current_parents: Set[str] = set()

    def flush():
        nonlocal current_id, current_name, current_parents
        if current_id and current_name:
            dag[current_id] = GOTerm(go_id=current_id, name=current_name, parents=set(current_parents))
        current_id, current_name, current_parents = None, None, set()

    with go_obo_path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line == "[Term]":
                flush()
                continue
            if not line or line.startswith("["):
                continue
            if line.startswith("id: GO:"):
                current_id = line.split("id:", 1)[1].strip()
            elif line.startswith("name:"):
                current_name = line.split("name:", 1)[1].strip()
            elif line.startswith("is_a:"):
                # is_a: GO:0008150 ! biological_process
                parent = line.split("is_a:", 1)[1].strip().split()[0]
                if parent.startswith("GO:"):
                    current_parents.add(parent)

    flush()

    # Compute depths (distance to root via is_a). Deterministic DP.
    def compute_depth(go_id: str, visiting: Set[str]) -> int:
        if go_id not in dag:
            return 0
        term = dag[go_id]
        if term.depth > 0:
            return term.depth
        if go_id in visiting:
            return 0
        visiting.add(go_id)
        if not term.parents:
            term.depth = 1
        else:
            parent_depths = [compute_depth(p, visiting) for p in term.parents if p in dag]
            term.depth = (max(parent_depths) + 1) if parent_depths else 1
        visiting.remove(go_id)
        return term.depth

    for gid in list(dag.keys()):
        compute_depth(gid, set())

    return dag


def load_goslim_ids(goslim_obo_path: str | Path) -> Set[str]:
    """Load GO-slim IDs from a goslim OBO file."""
    slim_dag = load_go_dag(goslim_obo_path)
    return set(slim_dag.keys())


def all_ancestors(go_id: str, dag: Dict[str, GOTerm]) -> Set[str]:
    """Return all is_a ancestors of go_id (excluding itself)."""
    if go_id not in dag:
        return set()
    out: Set[str] = set()
    stack = list(dag[go_id].parents)
    while stack:
        p = stack.pop()
        if p in out:
            continue
        out.add(p)
        if p in dag:
            stack.extend(list(dag[p].parents))
    return out


def project_to_goslim(go_id: str, dag: Dict[str, GOTerm], goslim_ids: Set[str]) -> Optional[str]:
    """Project GO term to a single GO-slim ancestor.

    Deterministic selection:
      - consider all slim ancestors
      - choose the deepest slim (most specific). ties broken by GO id.

    Returns selected slim GO id or None.
    """
    if go_id not in dag:
        return None
    anc = all_ancestors(go_id, dag)
    slim_anc = [a for a in anc if a in goslim_ids]
    if not slim_anc:
        return None
    slim_anc.sort(key=lambda sid: (-dag[sid].depth if sid in dag else 0, sid))
    return slim_anc[0]


# ----------------------------
# I/O: GOA GAF parsing
# ----------------------------


def load_goa_gaf(
    gaf_path: str | Path,
    universe_genes: Set[str],
    go_aspects: Sequence[str] = ("P",),
    evidence_allowlist: Optional[Set[str]] = None,
) -> Dict[str, TermSet]:
    """Load GOA GAF and build GO term -> gene set mapping.

    - universe_genes restricts associations to your panel/universe.
    - go_aspects: subset of {'P','F','C'}.
    - evidence_allowlist: if set, keep only these evidence codes.

    Returns term_id -> TermSet with source set to GO:BP/GO:MF/GO:CC.
    """
    gaf_path = Path(gaf_path)
    go_aspects = tuple(a.upper() for a in go_aspects)

    term2genes: Dict[str, Set[str]] = defaultdict(set)
    term2aspect: Dict[str, str] = {}

    with gaf_path.open("r", encoding="utf-8") as f:
        for line in f:
            if not line or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            # GAF 2.x fields
            # 1 DB, 2 DB Object ID, 3 DB Object Symbol, 4 Qualifier, 5 GO ID,
            # 6 DB:Reference, 7 Evidence Code, 8 With/From, 9 Aspect
            symbol = norm_gene(parts[2])
            go_id = parts[4].strip()
            ev = parts[6].strip().upper()
            aspect = parts[8].strip().upper()

            if aspect not in go_aspects:
                continue
            if evidence_allowlist is not None and ev not in evidence_allowlist:
                continue
            if symbol not in universe_genes:
                continue
            if not go_id.startswith("GO:"):
                continue

            term2genes[go_id].add(symbol)
            term2aspect[go_id] = aspect

    out: Dict[str, TermSet] = {}
    for go_id, genes in term2genes.items():
        aspect = term2aspect.get(go_id, "P")
        source = {"P": "GO:BP", "F": "GO:MF", "C": "GO:CC"}.get(aspect, "GO:BP")
        # term_name will be filled later from OBO if available
        out[go_id] = TermSet(term_id=go_id, term_name=go_id, source=source, genes=set(genes))
    return out


def apply_go_names(termsets: Dict[str, TermSet], go_dag: Dict[str, GOTerm]) -> Dict[str, TermSet]:
    """Fill term_name from GO DAG where possible."""
    out: Dict[str, TermSet] = {}
    for tid, ts in termsets.items():
        name = ts.term_name
        if tid in go_dag:
            name = go_dag[tid].name
        out[tid] = TermSet(term_id=tid, term_name=name, source=ts.source, genes=set(ts.genes))
    return out


# ----------------------------
# Enrichment: ORA (hypergeometric)
# ----------------------------


def ora_hypergeom(
    markers: Set[str],
    universe: Set[str],
    termsets: Dict[str, TermSet],
    min_overlap: int = 2,
    max_set_size: int = 500,
) -> List[EnrichHit]:
    """Hypergeometric ORA across all provided termsets."""
    M = len(universe)
    n = len(markers)

    hits: List[EnrichHit] = []
    pvals: List[float] = []

    for term_id, ts in termsets.items():
        genes = ts.genes & universe
        K = len(genes)
        if K < 2 or K > max_set_size:
            continue
        overlap_genes = genes & markers
        k = len(overlap_genes)
        if k < min_overlap:
            continue
        # P(X >= k)
        p = hypergeom.sf(k - 1, M, K, n)
        hits.append(
            EnrichHit(
                term_id=term_id,
                term_name=ts.term_name,
                source=ts.source,
                set_size=K,
                overlap=k,
                overlap_genes=set(overlap_genes),
                pval=float(p),
                fdr=math.nan,
            )
        )
        pvals.append(float(p))

    if not hits:
        return []

    _, qvals, _, _ = multipletests(pvals, method="fdr_bh")
    hits = [
        EnrichHit(
            term_id=h.term_id,
            term_name=h.term_name,
            source=h.source,
            set_size=h.set_size,
            overlap=h.overlap,
            overlap_genes=h.overlap_genes,
            pval=h.pval,
            fdr=float(q),
        )
        for h, q in zip(hits, qvals)
    ]
    hits.sort(key=lambda x: (x.fdr, x.pval, -x.overlap))
    return hits


# ----------------------------
# Redundancy reduction: term clustering
# ----------------------------


def jaccard(a: Set[str], b: Set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def cluster_terms(hits: List[EnrichHit], sim_threshold: float = 0.35, max_terms: int = 40) -> List[List[EnrichHit]]:
    """Greedy clustering of terms by gene-overlap Jaccard."""
    clusters: List[List[EnrichHit]] = []
    for h in hits[:max_terms]:
        placed = False
        for cl in clusters:
            rep = cl[0]
            if jaccard(h.overlap_genes, rep.overlap_genes) >= sim_threshold:
                cl.append(h)
                placed = True
                break
        if not placed:
            clusters.append([h])
    return clusters


def choose_representative(cluster: List[EnrichHit]) -> EnrichHit:
    """Choose the most significant term in a cluster (deterministic)."""
    cluster_sorted = sorted(cluster, key=lambda x: (x.fdr, x.pval, -x.overlap, x.term_id))
    return cluster_sorted[0]


# ----------------------------
# Program inference
# ----------------------------


def extract_go_id(term_id: str) -> Optional[str]:
    """Return GO id if term_id is or contains GO:#######."""
    if term_id.startswith("GO:") and _RE_GO_ID.fullmatch(term_id):
        return term_id
    m = _RE_GO_ID.search(term_id)
    return m.group(0) if m else None


def build_programs(
    clusters: List[List[EnrichHit]],
    go_dag: Optional[Dict[str, GOTerm]] = None,
    goslim_ids: Optional[Set[str]] = None,
    max_programs: int = 8,
) -> List[Dict]:
    """Build functional programs.

    If go_dag + goslim_ids are provided, GO terms are projected to GO-slim and
    programs are keyed by GO-slim term name.

    Otherwise, programs correspond to clustered representative terms.
    """

    if go_dag is not None and goslim_ids is not None:
        # GO-slim mode (tissue-agnostic, no hardcoded rules):
        #   - Project *each* GO term to a GO-slim ancestor
        #   - Group evidence terms by GO-slim label
        #   - Non-GO pathways (e.g., Reactome/Hallmark) are kept as-is (cluster representative label)
        program2genes: Dict[str, Set[str]] = defaultdict(set)
        program2terms: Dict[str, List[Dict]] = defaultdict(list)

        for cl in clusters:
            rep = choose_representative(cl)
            non_go_label = rep.term_name

            # For each enriched term in the cluster, attempt GO-slim projection.
            # If the term is not GO-resolvable, attach it to the cluster's non-GO label.
            for h in cl:
                go_id = extract_go_id(h.term_id)
                if go_id:
                    slim_id = project_to_goslim(go_id, go_dag, goslim_ids)
                    if slim_id and slim_id in go_dag:
                        label = go_dag[slim_id].name
                    else:
                        # If no slim ancestor, keep the original GO term name (still GO-grounded)
                        label = h.term_name
                else:
                    label = non_go_label

                program2genes[label].update(h.overlap_genes)
                program2terms[label].append(
                    {
                        "term_id": h.term_id,
                        "term_name": h.term_name,
                        "source": h.source,
                        "pval": float(h.pval),
                        "fdr": float(h.fdr),
                        "overlap": int(h.overlap),
                        "set_size": int(h.set_size),
                        "overlap_genes": sorted(h.overlap_genes),
                    }
                )

        def rk(label: str):
            fdrs = [t["fdr"] for t in program2terms[label]]
            best = min(fdrs) if fdrs else 1.0
            return (best, -len(program2genes[label]), label)

        out: List[Dict] = []
        for label in sorted(program2genes.keys(), key=rk)[:max_programs]:
            out.append(
                {
                    "function": label,
                    "supporting_genes": sorted(program2genes[label]),
                    "supporting_terms": sorted(program2terms[label], key=lambda d: (d["fdr"], -d["overlap"], d["term_name"])),
                }
            )
        return out

    # Fallback: clustered representative terms
    out2: List[Dict] = []
    for cl in clusters[:max_programs]:
        rep = choose_representative(cl)
        genes = set().union(*[h.overlap_genes for h in cl])
        out2.append(
            {
                "function": rep.term_name,
                "supporting_genes": sorted(genes),
                "supporting_terms": [
                    {
                        "term_id": h.term_id,
                        "term_name": h.term_name,
                        "source": h.source,
                        "pval": float(h.pval),
                        "fdr": float(h.fdr),
                        "overlap": int(h.overlap),
                        "set_size": int(h.set_size),
                        "overlap_genes": sorted(h.overlap_genes),
                    }
                    for h in sorted(cl, key=lambda x: (x.fdr, x.pval, -x.overlap, x.term_id))
                ],
            }
        )
    return out2


def functional_programs_for_marker_set(
    marker_genes: Sequence[str],
    universe_genes: Sequence[str],
    *,
    gmt_paths: Optional[Sequence[str | Path]] = None,
    use_goa: bool = False,
    goa_gaf: Optional[str | Path] = None,
    go_obo: Optional[str | Path] = None,
    go_aspects: Sequence[str] = ("P",),
    go_evidence: Optional[Set[str]] = None,
    fdr_cutoff: float = 0.10,
    min_overlap: int = 2,
    max_set_size: int = 500,
    sim_threshold: float = 0.35,
    max_programs: int = 8,
    goslim_obo: Optional[str | Path] = None,
) -> List[Dict]:
    markers = {norm_gene(g) for g in marker_genes if g and str(g).strip()}
    universe = {norm_gene(g) for g in universe_genes if g and str(g).strip()}

    go_dag: Optional[Dict[str, GOTerm]] = None
    goslim_ids: Optional[Set[str]] = None

    if go_obo is not None:
        go_dag = load_go_dag(go_obo)
    if goslim_obo is not None:
        if go_dag is None:
            raise ValueError("--goslim-obo requires --go-obo")
        goslim_ids = load_goslim_ids(goslim_obo)

    # Build termsets
    termsets: Dict[str, TermSet] = {}

    if use_goa:
        if goa_gaf is None or go_obo is None:
            raise ValueError("--use-goa requires --goa-gaf and --go-obo")
        termsets = load_goa_gaf(goa_gaf, universe_genes=universe, go_aspects=go_aspects, evidence_allowlist=go_evidence)
        if go_dag is not None:
            termsets = apply_go_names(termsets, go_dag)
    else:
        if not gmt_paths:
            raise ValueError("Provide --gmt paths or set --use-goa")
        for p in gmt_paths:
            termsets.update(load_gmt(p))

    hits = ora_hypergeom(markers=markers, universe=universe, termsets=termsets, min_overlap=min_overlap, max_set_size=max_set_size)
    hits = [h for h in hits if h.fdr <= fdr_cutoff]

    if not hits:
        return [
            {
                "function": "Unresolved functional enrichment",
                "supporting_genes": sorted(markers),
                "supporting_terms": [],
            }
        ]

    clusters = cluster_terms(hits, sim_threshold=sim_threshold, max_terms=40)
    return build_programs(clusters, go_dag=go_dag, goslim_ids=goslim_ids, max_programs=max_programs)


# ----------------------------
# CLI
# ----------------------------


def read_universe(universe_path: str | Path) -> List[str]:
    """Read universe genes from a one-column CSV/TSV."""
    p = Path(universe_path)
    df = pd.read_csv(p)
    if df.shape[1] < 1:
        raise ValueError(f"Universe file has no columns: {p}")
    return df.iloc[:, 0].astype(str).tolist()


def read_marker_table(markers_path: str | Path, cluster_col: str, gene_col: str, score_col: Optional[str], top_n: int) -> Dict[int, List[str]]:
    """Read marker genes and select top_n per cluster."""
    p = Path(markers_path)
    df = pd.read_csv(p)

    for col in [cluster_col, gene_col]:
        if col not in df.columns:
            raise ValueError(f"Markers file missing required column '{col}'. Columns: {list(df.columns)}")

    if score_col and score_col in df.columns:
        df = df.sort_values([cluster_col, score_col], ascending=[True, False])
    else:
        # stable deterministic order within each cluster
        df = df.sort_values([cluster_col, gene_col], ascending=[True, True])

    out: Dict[int, List[str]] = {}
    for cid, sub in df.groupby(cluster_col):
        genes = sub[gene_col].astype(str).head(top_n).tolist()
        try:
            cid_int = int(cid)
        except Exception:
            # if cluster ids are strings, preserve stable ordering by hashing is risky; just enumerate
            cid_int = int(float(cid)) if str(cid).replace('.', '', 1).isdigit() else None  # type: ignore
        out[int(cid_int) if cid_int is not None else cid] = genes  # type: ignore
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="Deterministic functional-program inference for clusters")

    ap.add_argument("--markers", required=True, help="CSV containing marker genes per cluster")
    ap.add_argument("--universe", required=True, help="CSV with one column listing universe/panel genes")
    ap.add_argument("--out", required=True, help="Output CSV path")

    ap.add_argument("--cluster-col", default="cluster", help="Cluster column name in markers CSV")
    ap.add_argument("--gene-col", default="gene", help="Gene column name in markers CSV")
    ap.add_argument("--score-col", default="log2FC", help="Score column for ranking markers (optional)")
    ap.add_argument("--top-n", type=int, default=20, help="Top N markers per cluster")

    # Enrichment source
    ap.add_argument("--gmt", action="append", default=[], help="GMT path (repeatable). Used when --use-goa is not set")

    ap.add_argument("--use-goa", action="store_true", help="Use GOA GAF + GO OBO backend for GO")
    ap.add_argument("--goa-gaf", default=None, help="GOA GAF file path (required with --use-goa)")
    ap.add_argument("--go-obo", default=None, help="GO OBO file path (recommended with --use-goa; required with --goslim-obo)")
    ap.add_argument("--go-aspects", default="P", help="GO aspects to include from GOA: comma-separated subset of P,F,C")
    ap.add_argument("--go-evidence", default=None, help="Comma-separated evidence codes to allow (optional)")

    # Programs
    ap.add_argument("--goslim-obo", default=None, help="GO-slim OBO file; if provided with --go-obo, GO terms are projected to GO-slim programs")

    # Thresholds
    ap.add_argument("--fdr", type=float, default=0.10, help="FDR cutoff")
    ap.add_argument("--min-overlap", type=int, default=2, help="Minimum overlap between markers and a term")
    ap.add_argument("--max-set-size", type=int, default=500, help="Maximum term gene-set size")
    ap.add_argument("--sim", type=float, default=0.35, help="Jaccard similarity threshold for term clustering")
    ap.add_argument("--max-programs", type=int, default=8, help="Max programs per cluster")

    args = ap.parse_args()

    universe_genes = read_universe(args.universe)

    score_col = args.score_col if args.score_col and args.score_col != "None" else None
    markers_by_cluster = read_marker_table(args.markers, args.cluster_col, args.gene_col, score_col, args.top_n)

    go_aspects = tuple(a.strip().upper() for a in str(args.go_aspects).split(",") if a.strip())
    go_evidence = None
    if args.go_evidence:
        go_evidence = {e.strip().upper() for e in str(args.go_evidence).split(",") if e.strip()}

    # Output
    out_rows: List[Dict[str, str]] = []

    # Ensure stable cluster ordering
    cluster_ids = sorted(markers_by_cluster.keys(), key=lambda x: (isinstance(x, str), x))

    for cid in cluster_ids:
        marker_genes = markers_by_cluster[cid]
        programs = functional_programs_for_marker_set(
            marker_genes,
            universe_genes,
            gmt_paths=args.gmt,
            use_goa=args.use_goa,
            goa_gaf=args.goa_gaf,
            go_obo=args.go_obo,
            go_aspects=go_aspects,
            go_evidence=go_evidence,
            fdr_cutoff=args.fdr,
            min_overlap=args.min_overlap,
            max_set_size=args.max_set_size,
            sim_threshold=args.sim,
            max_programs=args.max_programs,
            goslim_obo=args.goslim_obo,
        )

        for pr in programs:
            out_rows.append(
                {
                    "cluster": str(cid),
                    "function": pr["function"],
                    "supporting_genes": ";".join(pr["supporting_genes"]),
                    "supporting_terms_json": json.dumps(pr.get("supporting_terms", []), separators=(",", ":")),
                }
            )

        print(f"Finished cluster {cid}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(out_rows).to_csv(out_path, index=False)
    print(f"Saved to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
