#!/usr/bin/env python3
"""
00_validate_go_annotation.py - check that gene-level GO files use the same gene
IDs as the RNA-seq count tables (i.e. were made from the same gene models).

Run from this folder:  python3 00_validate_go_annotation.py
Standard library only. Paths are set in CONFIG below.

Checks, per species:
  1. GO-annotated gene IDs that are absent from the final annotation GFF (expect 0).
  2. % of genes with a GO annotation in each DESeq2 baseMean decile. With matching
     IDs this rises with expression (well-expressed genes are better annotated);
     with mismatched IDs it is flat.
  3. Median expression percentile of genes annotated to "cytosolic large ribosomal
     subunit" (GO:0022625): high (~85-90) when IDs match, ~50 when they do not.
  4. Median GO-set Jaccard similarity between 1:1 orthologs (OrthoFinder) of the
     reference species and each other species, versus randomly shuffled pairs.
     ~1.0 when IDs match, ~0.08 (= shuffled) when they do not.

Writes results/annotation_validation.tsv.
"""

import csv
import gzip
import random
import re
import statistics
from pathlib import Path

GO_DIR = Path("go_files_excon_v2.5.0")
CONFIG = {
    # species: (GO file, DESeq2 results, final annotation GFF, OrthoFinder column)
    "BJ": (GO_DIR / "Belonogaster_juncea.go.txt", "../../res_Bjun 1.tsv",
           "../../../Input_genomes_and_gff/Belonogaster_juncea_final_annot.gff.gz",
           "Belonogaster_juncea"),
    "LF": (GO_DIR / "Liostenogaster_flavolineata.go.txt", "../../res_LF.tsv",
           "../../../Input_genomes_and_gff/Liostenogaster_flavolineata_final_annot.gff.gz",
           "Liostenogaster_flavolineata"),
    "PD": (GO_DIR / "Polistes_dominula.go.txt", "../../res_Pdom.tsv",
           "../../../Input_genomes_and_gff/Polistes_dominula_final_annot.gff.gz",
           "Polistes_dominula"),
}
ORTHOGROUPS = "../../Plot_overlap/Orthogroups.tsv"
REFERENCE = "BJ"            # species whose GO file the others are compared against
RIBOSOME_GO = "GO:0022625"  # cytosolic large ribosomal subunit
SEED = 1


def open_any(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_go(path):
    gene2go = {}
    with open_any(path) as fh:
        for line in fh:
            gene, go = line.rstrip("\n").split("\t")[:2]
            gene2go.setdefault(gene, set()).add(go)
    return gene2go


def read_gff_genes(path):
    genes = set()
    with open_any(path) as fh:
        for line in fh:
            cols = line.split("\t")
            if len(cols) > 8 and cols[2] == "gene":
                m = re.search(r"ID=([^;\s]+)", cols[8])
                if m:
                    genes.add(m.group(1))
    return genes


def read_basemean(path):
    with open(path) as fh:
        return {r["gene"]: float(r["baseMean"]) for r in csv.DictReader(fh, delimiter="\t")}


def one_to_one(rows, col_a, col_b):
    """1:1 ortholog pairs (gene IDs, isoform suffix removed) between two OrthoFinder columns."""
    pairs = []
    for r in rows[1:]:
        a = [x.strip() for x in r[col_a].split(",") if x.strip()] if len(r) > col_a else []
        b = [x.strip() for x in r[col_b].split(",") if x.strip()] if len(r) > col_b else []
        if len(a) == 1 and len(b) == 1:
            pairs.append((a[0].split(".")[0], b[0].split(".")[0]))
    return pairs


def jaccard(x, y):
    return len(x & y) / len(x | y)


def main():
    random.seed(SEED)
    go = {sp: read_go(cfg[0]) for sp, cfg in CONFIG.items()}
    with open(ORTHOGROUPS) as fh:
        og_rows = list(csv.reader(fh, delimiter="\t"))
    og_col = {h: i for i, h in enumerate(og_rows[0])}

    out_rows = []
    for sp, (go_path, res_path, gff_path, og_name) in CONFIG.items():
        gene2go = go[sp]
        genes = read_gff_genes(gff_path)
        bm = read_basemean(res_path)

        ranked = sorted(bm, key=bm.get)
        n = len(ranked)
        deciles = []
        for d in range(10):
            chunk = ranked[d * n // 10:(d + 1) * n // 10]
            deciles.append(100 * sum(g in gene2go for g in chunk) / len(chunk))
        pct = {g: 100 * i / (n - 1) for i, g in enumerate(ranked)}
        rib = [pct[g] for g, s in gene2go.items() if RIBOSOME_GO in s and g in pct]

        jac_orth = jac_shuf = n_pairs = None
        if sp != REFERENCE:
            ref = go[REFERENCE]
            pairs = [(a, b) for a, b in one_to_one(og_rows, og_col[CONFIG[REFERENCE][3]], og_col[og_name])
                     if a in ref and b in gene2go]
            shuffled = [b for _, b in pairs]
            random.shuffle(shuffled)
            jac_orth = statistics.median(jaccard(ref[a], gene2go[b]) for a, b in pairs)
            jac_shuf = statistics.median(jaccard(ref[a], gene2go[b]) for (a, _), b in zip(pairs, shuffled))
            n_pairs = len(pairs)

        row = {
            "species": sp,
            "go_file": str(go_path),
            "annotated_genes": len(gene2go),
            "annotated_not_in_final_gff": len(set(gene2go) - genes),
            "annotated_in_count_table": len(set(gene2go) & set(bm)),
            "pct_annotated_by_expr_decile": " ".join(f"{x:.0f}" for x in deciles),
            "ribosome_LSU_median_expr_pct": round(statistics.median(rib), 1) if rib else "",
            "ribosome_LSU_genes": len(rib),
            f"orth_pairs_vs_{REFERENCE}": n_pairs if n_pairs is not None else "",
            "orth_GO_jaccard_median": round(jac_orth, 2) if jac_orth is not None else "",
            "shuffled_GO_jaccard_median": round(jac_shuf, 2) if jac_shuf is not None else "",
        }
        out_rows.append(row)
        print(f"{sp}: not in GFF={row['annotated_not_in_final_gff']}, "
              f"deciles={row['pct_annotated_by_expr_decile']}, "
              f"ribosome pct={row['ribosome_LSU_median_expr_pct']}, "
              f"ortholog Jaccard={row['orth_GO_jaccard_median']} (shuffled {row['shuffled_GO_jaccard_median']})")

    Path("results").mkdir(exist_ok=True)
    with open("results/annotation_validation.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print("Saved: results/annotation_validation.tsv")


if __name__ == "__main__":
    main()
