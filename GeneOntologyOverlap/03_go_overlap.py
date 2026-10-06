#!/usr/bin/env python3
"""
go_overlap.py - compare significant GO terms across caste-biased gene lists.

Reads TopGO-style enrichment tables (one per list), keeps the terms that pass a
chosen p-value column/threshold, and reports how many terms are shared between
lists, which ones, and an UpSet plot of the intersections. With --manifest it
also tests whether each pairwise overlap is larger than expected by chance.

Input (either):
  * one Excel workbook, one sheet per list (sheet name = list name), e.g.
        "B. juncea queen-biased", "P. dominula worker-biased", ...
    Sheets without a GO ID column (e.g. an "About" sheet) are skipped.
  * several CSV/TSV files (file name without extension = list name)

Expected columns (case-insensitive): GO ID (GO.ID or GO_ID), term name (Term or
GO_term), ontology, and the p-value column you choose. Column aliases:
  none (uncorrected P) : none, P_value, classicFisher
  BH                   : BH, FDR_BH
TopGO strings such as "< 1e-30" are handled. Terms are matched on GO ID, not on
the (possibly truncated) term name.

Species and caste are parsed from the list name: the species is everything
before the caste word; "worker" or "non(-)rep" = non-reproductive list,
"queen" or "rep" = reproductive list.

Overlap significance (--manifest results/overlap_manifest.tsv, written by
01_run_topGO.R; gives each list's DEG file and its species' background genes and
gene-to-GO-term membership exactly as topGO used them):
  * Permutation test (primary). GO terms are not independent: a gene annotated
    to a term is also annotated to all of its parent terms, so enriched terms
    come in correlated clusters. To keep that structure, each permutation draws
    random DEG sets of the observed sizes from each species' background (queen
    and worker sets drawn without overlap, as the real lists are), recomputes
    the topGO classic Fisher P-value of every term, and takes the same number of
    top-ranked terms as the observed list has (ties broken at random). The null
    distribution of the pairwise overlap comes from these random term sets, so
    it keeps the GO graph structure, term-size effects and the list sizes.
    P = (1 + #null overlaps >= observed) / (1 + n_perm).
  * Hypergeometric test (for reference; treats terms as independent, so it is
    anti-conservative). Universe = terms tested in both lists that could reach
    the threshold in both (minimum attainable Fisher P < alpha).
  Only pairs of lists from different species are tested (within a species the
  queen- and worker-biased lists are the two directions of one comparison).
  P-values of both tests are BH-adjusted across all tested pairs (*_BH; used in
  the plot) and, separately, across only the same-caste pairs (*_BH_same_caste),
  i.e. the convergence hypotheses, with the opposite-caste pairs as controls.

Usage:
  python 03_go_overlap.py Supplementary_Table_3_GO_enrichment.xlsx -o go_overlap --alpha 0.01 \\
      --manifest results/overlap_manifest.tsv
  python 03_go_overlap.py Supplementary_Table_3_GO_enrichment.xlsx --pcol BH --alpha 0.05
  python 03_go_overlap.py lists/*.tsv -o go_overlap

Outputs (in the output folder):
  summary.tsv                 per list: terms in input, significant (all/BP/MF/CC),
                              and how many are unique to that list
  pairwise_overlap.tsv        every pair of lists: comparison type, sizes,
                              shared terms (count, Jaccard, names); with
                              --manifest also the permutation and hypergeometric
                              test results
  pairwise_matrix.tsv         square matrix of shared-term counts
  all_significant_terms.tsv   every term significant in >= 1 list, with presence
                              and p-value in each list (supplementary-table ready)
  shared_terms.tsv            the subset significant in >= 2 lists
  threshold_sensitivity.tsv   how the overlap changes under other thresholds
                              (bonferroni / BH / uncorrected; columns absent
                              from the input are skipped)
  upset.pdf, upset.png        UpSet plot of the intersections
  overlap_stats.pdf, .png     (--manifest) observed pairwise overlap vs the
                              permutation null, with BH-adjusted significance
  overlap_figure.pdf, .png    (--manifest) both in one figure: (A) UpSet plot,
                              (B) overlap vs permutation null
  permutation_check.tsv       (--manifest) recomputed vs topGO P-values for the
                              observed lists (should agree to rounding error)
  overlap_null_counts.tsv.gz  (--manifest) the permutation null: shared-term count
                              for every pair in every permutation

Requires: python3, pandas, openpyxl (for .xlsx), matplotlib; numpy and scipy
for --manifest
"""

import argparse
import itertools
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Colours follow the manuscript's Figure 6 convention (blue = reproductive,
# red = non-reproductive); change here if needed.
CASTE_COLOURS = {"rep": "#2f6db5", "non-rep": "#c8423b", "": "#7a7a7a"}
INK = "#2b2b2b"
MUTED = "#bdbdbd"
FAINT = "#f2f2f2"

# Extra thresholds reported in threshold_sensitivity.tsv (column, alpha).
SENSITIVITY = [("bonferroni", 0.05), ("BH", 0.05), ("none", 0.05), ("none", 0.01)]

# Accepted input column names for each p-value column.
P_ALIASES = {
    "none": ["none", "P_value", "classicFisher"],
    "BH": ["BH", "FDR_BH"],
}

# Permutations are run in chunks of this size (memory scales with the chunk).
PERM_CHUNK = 1000

# Order of comparison types in tables and the overlap_stats plot.
COMPARISON_ORDER = ["different species, same caste",
                    "different species, opposite caste",
                    "same species, opposite caste",
                    "same species"]


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #
def read_lists(paths):
    """Return {list_name: DataFrame} in input order."""
    lists = {}
    for p in map(Path, paths):
        if p.suffix.lower() in (".xlsx", ".xlsm", ".xls"):
            sheets = pd.read_excel(p, sheet_name=None)
            for name, df in sheets.items():
                if df.dropna(how="all").empty:
                    continue
                lists[name.strip()] = df
        else:
            sep = "\t" if p.suffix.lower() in (".tsv", ".txt") else ","
            lists[p.stem] = pd.read_csv(p, sep=sep)
    if not lists:
        sys.exit("No lists found in the input.")
    return lists


def find_col(df, *wanted):
    """Find a column by any of several names, ignoring case, dots, underscores and spaces."""
    norm = lambda s: re.sub(r"[\s._]", "", str(s)).lower()
    for w in wanted:
        for c in df.columns:
            if norm(c) == norm(w):
                return c
    return None


def to_p(series):
    """Numeric p-values; TopGO strings like '< 1e-30' become 1e-30."""
    s = series.astype(str).str.replace("<", "", regex=False).str.strip()
    return pd.to_numeric(s, errors="coerce")


def standardise(name, df, pcols):
    """Return a tidy frame: GO_ID, Term, ontology, and requested p columns.
    Returns None for sheets that are not enrichment tables (no GO ID column)."""
    id_col = find_col(df, "GO.ID", "GOID", "GO")
    if id_col is None:
        print(f"  skipping [{name}]: no GO ID column", file=sys.stderr)
        return None
    out = pd.DataFrame({"GO_ID": df[id_col].astype(str).str.strip()})
    term_col = find_col(df, "Term", "GO_term")
    ont_col = find_col(df, "ontology")
    out["Term"] = df[term_col].astype(str) if term_col else ""
    out["ontology"] = df[ont_col].astype(str).str.upper() if ont_col else ""
    for pc in pcols:
        c = find_col(df, *P_ALIASES.get(pc, [pc]))
        out[pc] = to_p(df[c]) if c is not None else float("nan")
    out = out[out["GO_ID"].str.startswith("GO:")]
    return out.drop_duplicates("GO_ID").reset_index(drop=True)


def parse_label(name):
    low = name.lower()
    if re.search(r"non[\s_-]*rep", low) or "nonrep" in low or "worker" in low:
        caste = "non-rep"
    elif "rep" in low or "queen" in low:
        caste = "rep"
    else:
        caste = ""
    species = re.split(r"[\s_]+(?:queen|worker|non|rep)", name, maxsplit=1,
                       flags=re.IGNORECASE)[0].strip()
    return species or name, caste


# --------------------------------------------------------------------------- #
# Overlap
# --------------------------------------------------------------------------- #
def significant_sets(tables, pcol, alpha):
    sets = {}
    for name, t in tables.items():
        if t[pcol].isna().all():
            print(f"  warning: [{name}] has no usable '{pcol}' values", file=sys.stderr)
        sets[name] = set(t.loc[t[pcol] < alpha, "GO_ID"])
    return sets


def comparison_type(a, b):
    sa, ca = parse_label(a)
    sb, cb = parse_label(b)
    if sa == sb:
        return "same species, opposite caste" if ca != cb else "same species"
    if ca == cb and ca:
        return "different species, same caste"
    return "different species, opposite caste"


def pairwise(sets, term_names):
    rows = []
    for a, b in itertools.combinations(sets, 2):
        shared = sets[a] & sets[b]
        union = sets[a] | sets[b]
        rows.append({
            "list_a": a,
            "list_b": b,
            "comparison": comparison_type(a, b),
            "n_a": len(sets[a]),
            "n_b": len(sets[b]),
            "shared": len(shared),
            "jaccard": round(len(shared) / len(union), 4) if union else 0.0,
            "shared_terms": "; ".join(
                f"{g} {term_names.get(g, '')}".strip() for g in sorted(shared)
            ),
        })
    return pd.DataFrame(rows)


def matrix(sets):
    names = list(sets)
    m = pd.DataFrame(0, index=names, columns=names)
    for a in names:
        for b in names:
            m.loc[a, b] = len(sets[a] & sets[b])
    return m


def term_table(tables, sets, pcol):
    """One row per term significant anywhere: presence + p-value per list."""
    all_terms = set().union(*sets.values())
    info = {}
    for t in tables.values():
        for r in t.itertuples(index=False):
            if r.GO_ID in all_terms and r.GO_ID not in info:
                info[r.GO_ID] = (r.Term, r.ontology)
    rows = []
    for g in all_terms:
        term, ont = info.get(g, ("", ""))
        row = {"GO_ID": g, "Term": term, "ontology": ont}
        in_lists = [n for n in sets if g in sets[n]]
        row["n_lists"] = len(in_lists)
        row["n_species"] = len({parse_label(n)[0] for n in in_lists})
        row["lists"] = "; ".join(in_lists)
        for n, t in tables.items():
            p = t.loc[t["GO_ID"] == g, pcol]
            row[f"{n} ({pcol})"] = p.iloc[0] if len(p) and g in sets[n] else None
        rows.append(row)
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["n_lists", "n_species", "GO_ID"],
                          ascending=[False, False, True]).reset_index(drop=True)


def summary(tables, sets):
    rows = []
    for n in sets:
        sp, caste = parse_label(n)
        others = set().union(*(sets[o] for o in sets if o != n))
        ont = tables[n].set_index("GO_ID")["ontology"]
        sig_ont = ont.reindex(sorted(sets[n]))
        rows.append({
            "list": n, "species": sp, "caste": caste,
            "terms_in_input": len(tables[n]),
            "significant": len(sets[n]),
            "BP": int((sig_ont == "BP").sum()),
            "MF": int((sig_ont == "MF").sum()),
            "CC": int((sig_ont == "CC").sum()),
            "unique_to_list": len(sets[n] - others),
        })
    return pd.DataFrame(rows)


def sensitivity(tables, thresholds):
    rows = []
    for pc, alpha in thresholds:
        if all(t[pc].isna().all() for t in tables.values()):
            continue
        sets = significant_sets(tables, pc, alpha)
        union = set().union(*sets.values())
        count = {g: sum(g in s for s in sets.values()) for g in union}
        species_of = {g: {parse_label(n)[0] for n in sets if g in sets[n]} for g in union}
        # distinct terms significant for the same caste in >= 2 species
        same_caste = set()
        for caste in ("rep", "non-rep"):
            for g in union:
                sp = {parse_label(n)[0] for n in sets
                      if parse_label(n)[1] == caste and g in sets[n]}
                if len(sp) >= 2:
                    same_caste.add(g)
        rows.append({
            "p_column": pc, "alpha": alpha,
            "terms_significant_in_any_list": len(union),
            "in_>=2_lists": sum(c >= 2 for c in count.values()),
            "in_>=2_species": sum(len(s) >= 2 for s in species_of.values()),
            "same_caste_in_>=2_species": len(same_caste),
            **{f"n_{n}": len(s) for n, s in sets.items()},
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Overlap significance (--manifest)
# --------------------------------------------------------------------------- #
def read_manifest(path, list_names):
    """{list_name: file paths} from 01_run_topGO.R's manifest (paths relative to it)."""
    m = pd.read_csv(path, sep="\t")
    base = Path(path).parent
    info = {r.list: {"species": r.species,
                     "results": base / r.results_file,
                     "degs": base / r.degs_file,
                     "membership": base / r.membership_file,
                     "background": base / r.background_file}
            for r in m.itertuples(index=False)}
    missing = [n for n in list_names if n not in info]
    if missing:
        sys.exit(f"Lists not found in manifest {path}: {missing}")
    return {n: info[n] for n in list_names}


class SpeciesGO:
    """Background genes and, per ontology, the gene x term incidence matrix as topGO used it."""

    def __init__(self, membership, background):
        from scipy import sparse
        self.genes = [l.strip() for l in open(background) if l.strip()]
        self.gene_idx = {g: i for i, g in enumerate(self.genes)}
        mem = pd.read_csv(membership, sep="\t")
        self.ont = {}
        for ont, df in mem.groupby("Ontology"):
            terms = pd.unique(df["GO_ID"])
            term_idx = {t: j for j, t in enumerate(terms)}
            rows = df["gene"].map(self.gene_idx)
            if rows.isna().any():
                sys.exit(f"{membership}: genes missing from background {background}")
            M = sparse.csr_matrix(
                (np.ones(len(df), dtype=np.int32),
                 (rows.to_numpy(dtype=np.int64), df["GO_ID"].map(term_idx).to_numpy())),
                shape=(len(self.genes), len(terms)))
            feasible = np.asarray(M.sum(axis=1)).ravel() > 0   # genes annotated in this ontology
            self.ont[ont] = {"terms": np.asarray(terms), "M": M,
                             "K": np.asarray(M.sum(axis=0)).ravel(),
                             "feasible": feasible, "N": int(feasible.sum())}

    def indicator(self, gene_sets):
        """Sparse (n_sets x n_genes) 0/1 matrix from lists of gene indices."""
        from scipy import sparse
        rows = np.repeat(np.arange(len(gene_sets)), [len(s) for s in gene_sets])
        cols = np.concatenate(gene_sets) if len(gene_sets) else np.array([], dtype=int)
        return sparse.csr_matrix((np.ones(len(cols), dtype=np.int32), (rows, cols)),
                                 shape=(len(gene_sets), len(self.genes)))

    def pvalues(self, X):
        """Classic Fisher (hypergeometric upper-tail) P for every term, one row per DEG set."""
        ps, terms = [], []
        for ont in sorted(self.ont):
            o = self.ont[ont]
            counts = (X @ o["M"]).toarray()
            n = np.asarray(X[:, np.flatnonzero(o["feasible"])].sum(axis=1)).ravel()
            ps.append(fisher_upper(counts, o["K"], o["N"], n))
            terms.append(o["terms"])
        return np.hstack(ps), np.concatenate(terms)

    def min_attainable_p(self, n_by_ont):
        """Smallest Fisher P each term could reach given the observed DEG count per ontology."""
        from scipy.stats import hypergeom
        out = {}
        for ont, o in self.ont.items():
            n = n_by_ont[ont]
            x = np.minimum(o["K"], n)
            out.update(zip(o["terms"], hypergeom.sf(x - 1, o["N"], o["K"], n)))
        return out


def fisher_upper(counts, K, N, n):
    """P(X >= counts) for X ~ Hypergeom(N, K, n): topGO's one-sided classic Fisher test.
    counts (R, T), K (T,), N scalar, n (R,). Each distinct (n, K, count) is computed once."""
    from scipy.stats import hypergeom
    R, T = counts.shape
    c_base, k_base = int(counts.max()) + 1, int(K.max()) + 1
    key = ((np.asarray(n, dtype=np.int64).reshape(-1, 1) * k_base + K.reshape(1, -1)) * c_base
           + counts).ravel()
    uniq, inv = np.unique(key, return_inverse=True)
    uc = uniq % c_base
    uk = (uniq // c_base) % k_base
    un = uniq // (c_base * k_base)
    return hypergeom.sf(uc - 1, N, uk, un)[inv].reshape(R, T)


def bh(p):
    """Benjamini-Hochberg adjusted P-values."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adj = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(adj)
    out[order] = np.minimum(adj, 1)
    return out


def overlap_tests(sets, manifest, alpha, n_perm, seed, outdir):
    """Permutation and hypergeometric tests for every pair of lists."""
    from scipy.stats import hypergeom
    rng = np.random.default_rng(seed)
    names = list(sets)
    species = {}
    for n in names:
        sp = manifest[n]["species"]
        if sp not in species:
            species[sp] = SpeciesGO(manifest[n]["membership"], manifest[n]["background"])

    # Observed DEG sets (background genes only), recomputed P-values, attainable P
    deg_idx, attainable, check = {}, {}, []
    for n in names:
        S = species[manifest[n]["species"]]
        degs = [l.strip() for l in open(manifest[n]["degs"]) if l.strip()]
        deg_idx[n] = np.array(sorted(S.gene_idx[g] for g in degs if g in S.gene_idx), dtype=np.int64)
        X = S.indicator([deg_idx[n]])
        p_obs, terms = S.pvalues(X)
        topgo = pd.read_csv(manifest[n]["results"], sep="\t", quoting=3,   # 3 = QUOTE_NONE
                            usecols=["GO_ID", "P_value"]).set_index("GO_ID")["P_value"]
        diff = np.abs(np.log10(p_obs[0]) - np.log10(topgo.reindex(terms).to_numpy()))
        check.append({"list": n, "terms": len(terms), "DEGs_in_background": len(deg_idx[n]),
                      "max_abs_log10P_diff_vs_topGO": float(np.nanmax(diff)),
                      "terms_missing_in_topGO": int(np.isnan(diff).sum())})
        n_by_ont = {ont: int(o["feasible"][deg_idx[n]].sum()) for ont, o in S.ont.items()}
        p_min = S.min_attainable_p(n_by_ont)
        attainable[n] = {t for t, p in p_min.items() if p < alpha}
        stray = sets[n] - set(terms)
        if stray:
            sys.exit(f"[{n}] {len(stray)} significant terms are not in the topGO graph")
    pd.DataFrame(check).to_csv(outdir / "permutation_check.tsv", sep="\t", index=False)

    # Permutation null: random top-k term sets per list (k = observed number significant),
    # in chunks of PERM_CHUNK permutations so memory stays flat
    all_terms = sorted(set().union(*(set(np.concatenate([o["terms"] for o in S.ont.values()]))
                                     for S in species.values())))
    gid = {t: i for i, t in enumerate(all_terms)}
    # Only between-species pairs are tested: within a species the queen- and worker-biased
    # lists are the two directions of one comparison, so their overlap is not a question.
    pairs = [(a, b) for a, b in itertools.combinations(names, 2)
             if manifest[a]["species"] != manifest[b]["species"]]
    sp_lists = {sp: [n for n in names if manifest[n]["species"] == sp] for sp in species}
    null_counts = {pr: [] for pr in pairs}
    done = 0
    while done < n_perm:
        size = min(PERM_CHUNK, n_perm - done)
        B = {}
        for sp, S in species.items():
            draws = {n: [] for n in sp_lists[sp]}
            for _ in range(size):
                perm = rng.permutation(len(S.genes))
                start = 0
                for n in sp_lists[sp]:                 # disjoint draws within a species
                    draws[n].append(perm[start:start + len(deg_idx[n])])
                    start += len(deg_idx[n])
            for n in sp_lists[sp]:
                k = len(sets[n])
                Bn = np.zeros((size, len(all_terms)), dtype=bool)
                if k:
                    P, terms = S.pvalues(S.indicator(draws[n]))
                    order = np.lexsort((rng.random(P.shape), P), axis=-1)[:, :k]
                    tcol = np.array([gid[t] for t in terms])
                    Bn[np.arange(size)[:, None], tcol[order]] = True
                B[n] = Bn
        for a, b in pairs:
            null_counts[(a, b)].append((B[a] & B[b]).sum(axis=1))
        done += size
        print(f"  permutations: {done}/{n_perm}", file=sys.stderr)
    if n_perm > 0:
        null_counts = {pr: np.concatenate(v) for pr, v in null_counts.items()}
        pd.DataFrame({f"{a} | {b}": v for (a, b), v in null_counts.items()}).to_csv(
            outdir / "overlap_null_counts.tsv.gz", sep="\t", index_label="permutation")

    rows = []
    for a, b in pairs:
        obs = len(sets[a] & sets[b])
        row = {"list_a": a, "list_b": b}
        if n_perm > 0:
            null = null_counts[(a, b)]
            row.update({
                "perm_null_mean": round(float(null.mean()), 3),
                "perm_null_2.5pct": float(np.percentile(null, 2.5)),
                "perm_null_97.5pct": float(np.percentile(null, 97.5)),
                "perm_fold": round(obs / null.mean(), 3) if null.mean() > 0 else np.nan,
                "perm_p": (1 + int((null >= obs).sum())) / (1 + n_perm),
            })
        U = attainable[a] & attainable[b]
        na, nb = len(sets[a] & U), len(sets[b] & U)
        x = len(sets[a] & sets[b] & U)
        exp = na * nb / len(U) if U else 0.0
        row.update({
            "hyper_universe": len(U),
            "hyper_expected": round(exp, 3),
            "hyper_fold": round(x / exp, 3) if exp > 0 else np.nan,
            "hyper_p": float(hypergeom.sf(x - 1, len(U), na, nb)) if U else 1.0,
        })
        rows.append(row)
    res = pd.DataFrame(rows)
    # BH across all pairs, and across only the between-species same-caste pairs
    # (the convergence hypotheses; the other pairs act as controls)
    same = np.array([comparison_type(a, b) == "different species, same caste" for a, b in pairs])
    for col in (["perm_p"] if n_perm > 0 else []) + ["hyper_p"]:
        res[f"{col}_BH"] = bh(res[col])
        res[f"{col}_BH_same_caste"] = np.nan
        if same.any():
            res.loc[same, f"{col}_BH_same_caste"] = bh(res.loc[same, col])
    return res


def stars(p):
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"


def label_tex(name):
    """List name with the species in italics (matplotlib mathtext)."""
    species, _ = parse_label(name)
    rest = name[len(species):].strip().replace("-biased", "")
    if species != name and re.fullmatch(r"[A-Za-z. ]+", species):
        return r"$\it{" + species.replace(" ", r"\ ") + "}$ " + rest
    return name


def stats_layout(pw):
    """Rows of the overlap_stats plot: a heading above each comparison type, then its pairs."""
    df = pw.rename(columns={"perm_null_2.5pct": "null_lo", "perm_null_97.5pct": "null_hi"})
    df = df[df["perm_p"].notna()].copy()           # tested (between-species) pairs only
    df["order"] = df["comparison"].map({c: i for i, c in enumerate(COMPARISON_ORDER)})
    df = df.sort_values(["order"], kind="stable").reset_index(drop=True)
    items, y = [], 0
    for comp, grp in df.groupby("order", sort=True):
        items.append(("head", COMPARISON_ORDER[comp], y))
        y += 1
        for r in grp.itertuples(index=False):
            items.append(("pair", r, y))
            y += 1
    return df, items, y


def stats_size(pw):
    """Figure size (inches) of the standalone overlap_stats plot."""
    return 8.2, 0.27 * stats_layout(pw)[2] + 1.0


def draw_overlap_stats(fig, pw, title, n_perm, rect=None):
    """Observed shared terms per pair vs the permutation null (mean and central 95%).
    Draws into fig (a Figure or SubFigure); rect = [left, bottom, width, height] of the axes."""
    import matplotlib.pyplot as plt

    df, items, n_rows = stats_layout(pw)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 8,
                         "axes.edgecolor": MUTED, "axes.linewidth": 0.6})
    ax = fig.add_axes(rect) if rect else fig.subplots()
    xmax = max(df["shared"].max(), df["null_hi"].max()) * 1.08 + 1
    p_floor = 1 / (1 + n_perm)
    for kind, obj, yy in items:
        if kind == "head":
            ax.text(-0.01, yy, obj[0].upper() + obj[1:], transform=ax.get_yaxis_transform(),
                    ha="right", va="center", fontsize=8, fontweight="bold", color=INK)
            continue
        r = obj
        ca, cb = parse_label(r.list_a)[1], parse_label(r.list_b)[1]
        colour = CASTE_COLOURS[ca] if (ca == cb and ca) else INK
        ax.plot([r.null_lo, r.null_hi], [yy, yy],
                color=MUTED, lw=5, solid_capstyle="butt", zorder=1)
        ax.plot([r.perm_null_mean] * 2, [yy - 0.28, yy + 0.28], color="#7a7a7a", lw=1.2, zorder=2)
        ax.scatter([r.shared], [yy], s=30, color=colour, zorder=3, linewidths=0)
        p_txt = f"P ≤ {p_floor:.3g}" if r.perm_p <= p_floor else f"P = {r.perm_p:.3g}"
        fold = f"{r.perm_fold:.1f}×" if pd.notna(r.perm_fold) else "–"
        ax.text(1.01, yy, f"{r.shared} vs {r.perm_null_mean:.1f}  ({fold})  {p_txt}  {stars(r.perm_p_BH)}",
                transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=7, color=INK)
        ax.text(-0.01, yy, f"{label_tex(r.list_a)}  vs  {label_tex(r.list_b)}",
                transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=7, color=INK)
    for kind, obj, yy in items:
        if kind == "head":
            ax.axhline(yy - 0.5, color=FAINT, lw=6, zorder=0)
    ax.text(1.01, -0.75, "observed vs expected (fold), permutation P, BH",
            transform=ax.get_yaxis_transform(), ha="left", va="center",
            fontsize=7, fontstyle="italic", color="#7a7a7a")
    ax.set_ylim(n_rows - 0.5, -1.0)
    ax.set_xlim(-0.6, xmax)
    ax.set_yticks([])
    ax.set_xlabel("Shared GO terms")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", color=FAINT, lw=0.8, zorder=0)
    ax.set_title(title, loc="left", fontsize=9, color=INK)
    handles = [
        plt.Line2D([], [], color=INK, marker="o", lw=0, markersize=5, label="Observed"),
        plt.Line2D([], [], color=MUTED, lw=5,
                   label=f"Permutation null, central 95% ({n_perm:,} permutations)"),
        plt.Line2D([], [], color="#7a7a7a", lw=1.2, label="Permutation null, mean"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.06 - 1.2 / n_rows),
              ncol=3, frameon=False, fontsize=7)


def save_figure(fig, out_stem):
    import matplotlib.pyplot as plt
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def overlap_stats_plot(pw, out_stem, title, n_perm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=stats_size(pw))
    draw_overlap_stats(fig, pw, title, n_perm)
    save_figure(fig, out_stem)


# --------------------------------------------------------------------------- #
# UpSet plot (drawn directly with matplotlib; no extra dependency)
# --------------------------------------------------------------------------- #
def upset_layout(sets, max_bars):
    """Intersections to draw and panel sizes in inches; None if no term is significant."""
    names = list(sets)
    union = set().union(*sets.values())
    if not union:
        return None
    patterns = {}
    for g in union:
        key = tuple(g in sets[n] for n in names)
        patterns[key] = patterns.get(key, 0) + 1
    inter = sorted(patterns.items(), key=lambda kv: (-kv[1], -sum(kv[0])))[:max_bars]
    widths = [1.3,                                        # set-size bars
              0.075 * max(len(n) for n in names) + 0.2,  # list names
              max(0.32 * len(inter), 1.5)]               # matrix
    heights = [1.7, 0.3 * len(names)]                    # intersection bars, matrix
    return names, inter, widths, heights


def upset_size(sets, max_bars):
    """Figure size (inches) of the standalone UpSet plot."""
    _, _, widths, heights = upset_layout(sets, max_bars)
    return sum(widths), sum(heights) + 0.6


def draw_upset(fig, sets, title, max_bars=40):
    """UpSet plot drawn into fig (a Figure or SubFigure)."""
    import matplotlib.pyplot as plt

    names, inter, widths, heights = upset_layout(sets, max_bars)
    n_rows, n_cols = len(names), len(inter)
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 8,
                         "axes.edgecolor": MUTED, "axes.linewidth": 0.6})
    gs = fig.add_gridspec(2, 3, width_ratios=widths, height_ratios=heights,
                          wspace=0.0, hspace=0.05)
    ax_bar = fig.add_subplot(gs[0, 2])
    ax_mat = fig.add_subplot(gs[1, 2], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_mat)
    ax_lab = fig.add_subplot(gs[1, 1], sharey=ax_mat)
    ax_leg = fig.add_subplot(gs[0, 0:2])
    ax_leg.axis("off")

    xs = range(n_cols)
    counts = [c for _, c in inter]
    shared = [sum(k) > 1 for k, _ in inter]
    ax_bar.bar(xs, counts, width=0.62,
               color=[INK if s else MUTED for s in shared])
    for x, c in zip(xs, counts):
        ax_bar.text(x, c, str(c), ha="center", va="bottom", fontsize=7, color=INK)
    ax_bar.set_ylabel("Terms in intersection")
    ax_bar.spines[["top", "right", "bottom"]].set_visible(False)
    ax_bar.tick_params(axis="x", bottom=False, labelbottom=False)
    ax_bar.set_ylim(0, max(counts) * 1.15)
    if title:
        ax_bar.set_title(title, loc="left", fontsize=9, color=INK)

    # matrix: row r (top = first list)
    ys = list(range(n_rows))[::-1]
    for r, y in enumerate(ys):
        if r % 2 == 0:
            ax_mat.axhspan(y - 0.5, y + 0.5, color=FAINT, zorder=0)
    for x, (key, _) in zip(xs, inter):
        on = [ys[r] for r, v in enumerate(key) if v]
        ax_mat.scatter([x] * n_rows, ys, s=34, color=MUTED, zorder=2, linewidths=0)
        ax_mat.scatter([x] * len(on), on, s=40, color=INK, zorder=3, linewidths=0)
        if len(on) > 1:
            ax_mat.plot([x, x], [min(on), max(on)], color=INK, lw=1.6, zorder=2)
    ax_mat.set_xlim(-0.6, n_cols - 0.4)
    ax_mat.set_ylim(-0.5, n_rows - 0.5)
    ax_mat.axis("off")

    sizes = [len(sets[n]) for n in names]
    colours = [CASTE_COLOURS[parse_label(n)[1]] for n in names]
    ax_set.barh(ys, sizes, height=0.55, color=colours)
    ax_set.tick_params(axis="y", left=False, labelleft=False)
    ax_set.set_xlabel("Terms per list")
    ax_set.spines[["top", "left", "right"]].set_visible(False)
    for y, s in zip(ys, sizes):
        ax_set.text(s, y, f"{s} ", ha="right", va="center", fontsize=7, color=INK)
    ax_set.set_xlim(max(sizes) * 1.3 if max(sizes) else 1, 0)

    # list names in their own column, between the set-size bars and the matrix
    for r, y in enumerate(ys):
        if r % 2 == 0:
            ax_lab.axhspan(y - 0.5, y + 0.5, color=FAINT, zorder=0)
        species, _ = parse_label(names[r])
        rest = names[r][len(species):].strip()
        if species != names[r] and re.fullmatch(r"[A-Za-z. ]+", species):
            label = r"$\it{" + species.replace(" ", r"\ ") + "}$ " + rest   # italicise species
        else:
            label = names[r]
        ax_lab.text(0.04, y, label, ha="left", va="center", color=INK)
    ax_lab.set_xlim(0, 1)
    ax_lab.axis("off")

    if any(parse_label(n)[1] for n in names):
        handles = [plt.Rectangle((0, 0), 1, 1, color=CASTE_COLOURS[c])
                   for c in ("rep", "non-rep")]
        ax_leg.legend(handles, ["Reproductive-biased", "Non-reproductive-biased"],
                      loc="lower left", frameon=False, fontsize=7)


def upset_plot(sets, out_stem, title, max_bars=40):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if upset_layout(sets, max_bars) is None:
        print("  no significant terms - UpSet plot skipped", file=sys.stderr)
        return
    fig = plt.figure(figsize=upset_size(sets, max_bars))
    draw_upset(fig, sets, title, max_bars)
    save_figure(fig, out_stem)


def combined_figure(sets, pw, out_stem, upset_title, stats_title, n_perm, max_bars=40):
    """One figure: (A) UpSet plot, (B) pairwise overlap vs the permutation null."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if upset_layout(sets, max_bars) is None:
        return
    h_a = upset_size(sets, max_bars)[1]
    h_b = stats_size(pw)[1] + 0.4
    fig = plt.figure(figsize=(11.0, h_a + h_b))
    sub_a, sub_b = fig.subfigures(2, 1, height_ratios=[h_a, h_b])
    draw_upset(sub_a, sets, upset_title, max_bars)
    draw_overlap_stats(sub_b, pw, stats_title, n_perm, rect=[0.33, 0.16, 0.40, 0.76])
    for sub, letter in ((sub_a, "A"), (sub_b, "B")):
        sub.text(0.005, 0.99, letter, ha="left", va="top", fontsize=13,
                 fontweight="bold", color=INK)
    save_figure(fig, out_stem)


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="xlsx workbook or CSV/TSV files")
    ap.add_argument("-o", "--outdir", default="go_overlap")
    ap.add_argument("--pcol", default="none",
                    help="p-value column defining significance "
                         "(default: none = uncorrected; e.g. BH, bonferroni)")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--max-bars", type=int, default=40,
                    help="maximum intersections drawn in the UpSet plot")
    ap.add_argument("--manifest",
                    help="results/overlap_manifest.tsv from 01_run_topGO.R; enables the "
                         "overlap significance tests")
    ap.add_argument("--n-perm", type=int, default=10000,
                    help="permutations for the overlap test (0 = hypergeometric only)")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    raw = read_lists(args.inputs)
    pcols = sorted({args.pcol, *(c for c, _ in SENSITIVITY)})
    tables = {n: standardise(n, df, pcols) for n, df in raw.items()}
    tables = {n: t for n, t in tables.items() if t is not None}
    if not tables:
        sys.exit("No enrichment tables (with a GO ID column) found in the input.")

    sets = significant_sets(tables, args.pcol, args.alpha)
    term_names = {}
    for t in tables.values():
        for g, term in zip(t["GO_ID"], t["Term"]):
            term_names.setdefault(g, term)

    summ = summary(tables, sets)
    pw = pairwise(sets, term_names)
    terms = term_table(tables, sets, args.pcol)
    shared = terms[terms["n_lists"] >= 2] if not terms.empty else terms
    thresholds = [(args.pcol, args.alpha)] + [t for t in SENSITIVITY
                                              if t != (args.pcol, args.alpha)]
    sens = sensitivity(tables, thresholds)
    sig_label = (f"P < {args.alpha}" if args.pcol == "none"
                 else f"{args.pcol} P < {args.alpha}")

    if args.manifest:
        manifest = read_manifest(args.manifest, list(sets))
        print("Testing pairwise overlaps ...", file=sys.stderr)
        tests = overlap_tests(sets, manifest, args.alpha, args.n_perm, args.seed, out)
        pw = pw.merge(tests, on=["list_a", "list_b"], how="left")
        cols = [c for c in pw.columns if c != "shared_terms"] + ["shared_terms"]
        pw = pw[cols]
        if args.n_perm > 0:
            stats_title = f"Shared enriched GO terms ({sig_label}) vs permutation null"
            overlap_stats_plot(pw, out / "overlap_stats", stats_title, args.n_perm)
            combined_figure(sets, pw, out / "overlap_figure",
                            f"Enriched GO terms ({sig_label})", stats_title,
                            args.n_perm, args.max_bars)

    summ.to_csv(out / "summary.tsv", sep="\t", index=False)
    pw.to_csv(out / "pairwise_overlap.tsv", sep="\t", index=False)
    matrix(sets).to_csv(out / "pairwise_matrix.tsv", sep="\t")
    terms.to_csv(out / "all_significant_terms.tsv", sep="\t", index=False)
    shared.to_csv(out / "shared_terms.tsv", sep="\t", index=False)
    sens.to_csv(out / "threshold_sensitivity.tsv", sep="\t", index=False)
    upset_plot(sets, out / "upset", f"Enriched GO terms ({sig_label})", args.max_bars)

    # short console report
    union = set().union(*sets.values())
    print(f"\nSignificance: {args.pcol} < {args.alpha}")
    print(summ[["list", "significant", "unique_to_list"]].to_string(index=False))
    print(f"\nTerms significant in any list: {len(union)}")
    print(f"Terms significant in >= 2 lists: {len(shared)}")
    if not shared.empty:
        print(f"Terms significant in >= 2 species: {int((shared['n_species'] >= 2).sum())}")
        print("\nShared terms:")
        for r in shared.itertuples(index=False):
            print(f"  {r.GO_ID}  {r.Term[:45]:45s}  {r.lists}")
    if args.manifest:
        show = ["list_a", "list_b", "comparison", "shared"]
        if args.n_perm > 0:
            show += ["perm_null_mean", "perm_fold", "perm_p", "perm_p_BH", "perm_p_BH_same_caste"]
        show += ["hyper_expected", "hyper_p", "hyper_p_BH"]
        print("\nPairwise overlap tests:")
        print(pw[show].to_string(index=False))
    print(f"\nOutputs written to {out.resolve()}")


if __name__ == "__main__":
    main()
