"""Summarise Tandem Repeats Finder output for the final B. juncea assembly (chromosomes 1-39).

Input:  B_juncea_trf.txt  (TRF v4.09.1, '2 7 7 80 10 50 2000 -h -ngs -l 10', one run per scaffold)
Also uses replot_omnic_final.py (Hi-C contacts per position), the final GFF and the excon N0.tsv for the gene overlap.
All input files are expected in this folder (see README.md).
Outputs: trf_summary.tsv, satellite_296bp_consensus.fa, printed summary.
Run from this folder with the hic-straw environment: python trf_summary.py
"""
import re, gzip, bisect, collections
import numpy as np
import pandas as pd

CHROMS = [f"HiC_scaffold_{i}" for i in range(1, 40)]
SAT_PERIODS = [(285, 305), (580, 600), (870, 900)]      # ~296-bp unit and arrays reported as 2 or 3 units
LONG = 10_000                                             # "long array" threshold (bp)

# --- read TRF records
recs, sc = [], None
for line in open("B_juncea_trf.txt"):
    if line.startswith("@"):
        sc = line[1:].strip(); continue
    f = line.split()
    recs.append((sc, int(f[0]), int(f[1]), int(f[2]), float(f[3]), int(f[5]), f[13]))
trf = pd.DataFrame(recs, columns=["scaf", "start", "end", "period", "copies", "pmatch", "consensus"])
trf = trf[trf.scaf.isin(CHROMS)].copy()
trf["len"] = trf.end - trf.start + 1
is_sat = np.zeros(len(trf), bool)
for a, b in SAT_PERIODS:
    is_sat |= trf.period.between(a, b).to_numpy()
trf["class"] = np.where(is_sat, "296-bp satellite",
               np.where(trf.period < 10, "microsatellite (<10 bp)",
               np.where(trf.period < 100, "minisatellite (10-99 bp)", "other satellite (>=100 bp)")))

seqlen, name = {}, None
for line in open("B_junceaGenome.FINAL.fasta"):
    if line.startswith(">"):
        name = line[1:].split()[0]; seqlen[name] = 0
    else:
        seqlen[name] += len(line.strip())
CHR_TOTAL = sum(seqlen[c] for c in CHROMS)


def merge(df):
    """Merged intervals per scaffold -> dict scaf -> list of (start, end)."""
    out = {}
    for s, g in df.groupby("scaf"):
        iv = sorted(zip(g.start, g.end)); m = [list(iv[0])]
        for a, b in iv[1:]:
            if a <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        out[s] = m
    return out


def bp(m):
    return sum(e - s + 1 for v in m.values() for s, e in v)


rows = []
all_m = merge(trf)
rows.append(("all tandem repeats", bp(all_m)))
for cl in ["296-bp satellite", "other satellite (>=100 bp)", "minisatellite (10-99 bp)", "microsatellite (<10 bp)"]:
    rows.append((cl, bp(merge(trf[trf["class"] == cl]))))
sat_m = merge(trf[trf["class"] == "296-bp satellite"])
rows.append(("296-bp satellite, arrays >= 10 kb", bp(merge(trf[(trf["class"] == "296-bp satellite") & (trf.len >= LONG)]))))
rows.append(("any class, arrays >= 10 kb", bp(merge(trf[trf.len >= LONG]))))
summ = pd.DataFrame(rows, columns=["class", "bp"])
summ["Mb"] = (summ.bp / 1e6).round(1)
summ["pct_of_chr1-39"] = (100 * summ.bp / CHR_TOTAL).round(1)
print(f"chromosomes 1-39: {CHR_TOTAL/1e6:.1f} Mb\n"); print(summ.to_string(index=False))
summ.to_csv("trf_summary.tsv", sep="\t", index=False)

sat = trf[trf["class"] == "296-bp satellite"]
print(f"\n296-bp satellite: {len(sat):,} arrays; {sat.scaf.nunique()} of 39 chromosomes; "
      f"arrays >= 100 kb: {(sat.len >= 1e5).sum()}, >= 1 Mb: {(sat.len >= 1e6).sum()}; longest {sat.len.max()/1e6:.2f} Mb "
      f"({sat.loc[sat.len.idxmax(), 'scaf']})")
per_chr = (pd.Series({c: sum(e - s + 1 for s, e in sat_m.get(c, [])) for c in CHROMS}) / 1e6)
print("chromosomes with most satellite (Mb):", per_chr.sort_values(ascending=False).head(8).round(1).to_dict())

# --- one family or several? compare each long array's monomer to the consensus of the longest 296-bp array
K = 10
rc = str.maketrans("ACGT", "TGCA")


def kmers(s):
    return {s[i:i + K] for i in range(len(s) - K + 1)}


ref = sat[sat.period.between(285, 305)].sort_values("len", ascending=False).iloc[0]
ref_seq = ref.consensus
ref_k = kmers(ref_seq + ref_seq) | kmers((ref_seq + ref_seq).translate(rc)[::-1])
with open("satellite_296bp_consensus.fa", "w") as fh:
    fh.write(f">Bjuncea_sat296 consensus of longest array {ref.scaf}:{ref.start}-{ref.end} period {ref.period}\n{ref_seq}\n")
long_sat = sat[sat.len >= LONG]
share = np.array([len(kmers(c) & ref_k) / max(1, len(kmers(c))) for c in long_sat.consensus])
w = long_sat.len.to_numpy()
print(f"\nlong (>= 10 kb) 296-bp arrays: {len(long_sat)}; sharing >= 50% of 10-mers with the reference monomer: "
      f"{100*np.mean(share >= 0.5):.0f}% of arrays, {100*w[share >= 0.5].sum()/w.sum():.0f}% of their length")
gc = (ref_seq.count("G") + ref_seq.count("C")) / len(ref_seq)
print(f"reference monomer: {len(ref_seq)} bp, GC {100*gc:.0f}% (written to satellite_296bp_consensus.fa)")

# --- overlap with Hi-C coverage (50-kb windows) and genes
src = open("replot_omnic_final.py").read()
exec(src.split("# --- read the original map")[0])
import hicstraw
mzd = hicstraw.HiCFile(HIC).getMatrixZoomData("assembly", "assembly", "observed", "NONE", "BP", READ_RES)
cov_old = np.zeros(n_old)
for a in range(0, L0, 20_000_000):
    b_ = min(a + 20_000_000, L0); rec = mzd.getRecords(a, b_ - 1, a, L0 - 1)
    x = np.fromiter((r.binX for r in rec), np.int64, len(rec)) // READ_RES
    y = np.fromiter((r.binY for r in rec), np.int64, len(rec)) // READ_RES
    v = np.fromiter((r.counts for r in rec), float, len(rec)); del rec
    m = (x <= y) & (x >= a // READ_RES) & (x < int(np.ceil(b_ / READ_RES))); x, y, v = x[m], y[m], v[m]
    np.add.at(cov_old, x, v); np.add.at(cov_old, y[x != y], v[x != y])
rel = cov_old / np.median(cov_old[cov_old > 0])
seg_k = []
for cn, s, e in chrom_bounds:
    k = 0
    for sg in segments:
        if s <= sg[3] < e:
            seg_k.append((sg[3], cn, s, k)); k += 1
seg_k.sort(); starts_ = [t[0] for t in seg_k]
W = 50_000
win_cov = collections.defaultdict(list)
for i in np.where(new_pos >= 0)[0]:
    j = bisect.bisect_right(starts_, new_pos[i]) - 1; _, cn, cs, k = seg_k[j]
    win_cov[(cn, (new_pos[i] - cs + GAP * k) // W)].append(rel[i])
win_cls = {key: ("no contacts (<10%)" if np.median(v) < 0.1 else "low (10-50%)" if np.median(v) < 0.5 else "normal (>50%)")
           for key, v in win_cov.items()}
sat_in_win = collections.Counter()
for c, ivs in sat_m.items():
    for s, e in ivs:
        for wi in range(s // W, e // W + 1):
            ov = min(e, (wi + 1) * W) - max(s, wi * W)
            if ov > 0:
                sat_in_win[win_cls.get((c, wi), "unknown")] += ov
tot_sat = sum(sat_in_win.values())
print("\n296-bp satellite by Hi-C coverage of its 50-kb window: "
      + ", ".join(f"{k} {100*v/tot_sat:.0f}%" for k, v in sat_in_win.most_common()))
nowin = [key for key, cl in win_cls.items() if cl == "no contacts (<10%)"]
def frac_sat(c, wi):
    return sum(max(0, min(e, (wi + 1) * W) - max(s, wi * W)) for s, e in sat_m.get(c, [])) / W
fs = np.array([frac_sat(*k) for k in nowin])
print(f"no-contact windows ({len(nowin)}, {len(nowin)*W/1e6:.0f} Mb): mean share covered by the 296-bp satellite {100*fs.mean():.0f}%")

# genes overlapping satellite arrays
n3 = pd.read_csv("N0.tsv", sep="\t", dtype=str).set_index("Orthogroup")
gl = lambda v: [] if pd.isna(v) else [g.strip() for g in v.split(",") if g.strip()]
others = [c for c in n3.columns if c != "Belonogaster_juncea"]
amp = {g for og, r in n3.iterrows()
       if len(gl(r.Belonogaster_juncea)) >= 20 and len(gl(r.Belonogaster_juncea)) >= 5 * max(1, max(len(gl(r[c])) for c in others))
       for g in gl(r.Belonogaster_juncea)}
sat_starts = {c: [s for s, e in v] for c, v in sat_m.items()}
def within(c, pos, flank):
    v = sat_m.get(c, []); i = bisect.bisect_right(sat_starts.get(c, []), pos) - 1
    for j in (i, i + 1):
        if 0 <= j < len(v) and v[j][0] - flank <= pos <= v[j][1] + flank:
            return True
    return False
res = collections.defaultdict(lambda: [0, 0, 0])
for line in gzip.open("Belonogaster_juncea_final_annot.gff.gz", "rt"):
    f = line.split("\t")
    if len(f) == 9 and f[2] == "gene" and f[0] in CHROMS:
        g = re.search(r"ID=([^;\s]+)", f[8]).group(1); mid = (int(f[3]) + int(f[4])) // 2
        key = "amplified (TE-like) families" if g in amp else "all other genes"
        res[key][0] += 1; res[key][1] += within(f[0], mid, 0); res[key][2] += within(f[0], mid, 10_000)
print("\ngenes inside 296-bp satellite arrays (or within 10 kb of one):")
for k, (n, a, b) in res.items():
    print(f"  {k:30s} n={n:5d}  inside {100*a/n:.1f}%  within 10 kb {100*b/n:.1f}%")
