"""Omni-C contact map of the final B. juncea assembly (Figure 1b), chromosomes 1-39 only.

Inputs (from Rick Masonbrink's 3D-DNA / Juicebox curation):
  B_junceaGenome.0.hic              contact map, drawn in the ORIGINAL (.0.assembly) layout
  B_junceaGenome.0.assembly         that original layout
  B_junceaGenome.0.review.assembly  the curated layout after Juicebox review
  B_junceaGenome.FINAL.fasta        final assembly (HiC_scaffold_1..99)

Each 50-kb bin of the original map is moved to its position in the reviewed layout (handling
split fragments and inversions), the reviewed scaffolds are matched to the final HiC_scaffold
names by length (3D-DNA adds a 500-bp gap at each join), and the 39 chromosomes are drawn in
numerical order at 250 kb. Needs: hic-straw (pip install hic-straw), numpy, matplotlib.
All input files are expected in this folder (see README.md).

Run from this folder:  python replot_omnic_final.py
"""
import re
import numpy as np
import hicstraw
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HIC, ASM0, ASMR = "B_junceaGenome.0.hic", "B_junceaGenome.0.assembly", "B_junceaGenome.0.review.assembly"
FASTA = "B_junceaGenome.FINAL.fasta"
READ_RES, PLOT_RES, N_CHROM, GAP = 50_000, 250_000, 39, 500


def read_assembly(path):
    frag, order = {}, []
    for line in open(path):
        if line.startswith(">"):
            name, idx, length = line[1:].split()
            frag[int(idx)] = (name, int(length))
        elif line.strip():
            order.append([int(x) for x in line.split()])
    return frag, order


def fasta_lengths(path):
    lens, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]; lens[name] = 0
        else:
            lens[name] += len(line.strip())
    return lens


# --- original layout: where each original fragment sits in the .hic, and its orientation
f0, o0 = read_assembly(ASM0)
orig = {}                       # name -> (start in map, length, orientation +1/-1)
pos = 0
for scaf in o0:
    for i in scaf:
        name, length = f0[abs(i)]
        orig[name] = (pos, length, 1 if i > 0 else -1); pos += length

# --- reviewed layout: each piece = (original fragment, native offset, length)
fr, orr = read_assembly(ASMR)
piece = {}
by_parent = {}
for idx, (name, length) in fr.items():
    if name in orig:
        piece[idx] = (name, 0, length)
    else:
        parent = name.split(":::")[0]
        k = int(re.search(r":::fragment_(\d+)", name).group(1))
        by_parent.setdefault(parent, []).append((k, idx, length))
for parent, parts in by_parent.items():
    off = 0
    for k, idx, length in sorted(parts):
        piece[idx] = (parent, off, length); off += length
    assert off == orig[parent][1], parent

# --- match reviewed scaffolds to final HiC_scaffold names by length
final = fasta_lengths(FASTA)
chrom_names = [f"HiC_scaffold_{i}" for i in range(1, N_CHROM + 1)]
unused = dict(final)
rev_name = {}
for k, scaf in enumerate(orr):
    raw = sum(fr[abs(i)][1] for i in scaf)
    for L in (raw + GAP * (len(scaf) - 1), raw):
        hit = [n for n, l in unused.items() if l == L]
        if hit:
            rev_name[k] = hit[0]; unused.pop(hit[0]); break
for k, scaf in enumerate(orr):          # leftovers: nearest remaining length (gap count can differ)
    if k not in rev_name and unused:    # the layout has 100 scaffolds, the FASTA 99: one tiny piece stays unnamed
        raw = sum(fr[abs(i)][1] for i in scaf)
        n = min(unused, key=lambda n: abs(unused[n] - raw)); rev_name[k] = n; unused.pop(n)
name_to_rev = {v: k for k, v in rev_name.items()}

# --- new coordinates: chromosomes 1..39 concatenated (without gaps), pieces in reviewed order
segments = []                   # (parent, native_start, length, new_start, new_orientation)
chrom_bounds, new = [], 0
for cn in chrom_names:
    start = new
    for i in orr[name_to_rev[cn]]:
        parent, off, length = piece[abs(i)]
        segments.append((parent, off, length, new, 1 if i > 0 else -1)); new += length
    chrom_bounds.append((cn, start, new))
total_new = new

# --- for each 50-kb bin of the original map, its position in the new layout (or -1)
L0 = sum(l for _, l in f0.values())
n_old = int(np.ceil(L0 / READ_RES))
old_centre = np.arange(n_old) * READ_RES + READ_RES // 2
new_pos = np.full(n_old, -1, dtype=np.int64)
seg_by_parent = {}
for s in segments:
    seg_by_parent.setdefault(s[0], []).append(s)
for parent, (ostart, olen, osign) in orig.items():
    if parent not in seg_by_parent:
        continue
    sel = np.where((old_centre >= ostart) & (old_centre < ostart + olen))[0]
    x = old_centre[sel] - ostart
    native = x if osign > 0 else olen - 1 - x
    for _, noff, nlen, nstart, nsign in seg_by_parent[parent]:
        m = (native >= noff) & (native < noff + nlen)
        rel = native[m] - noff
        new_pos[sel[m]] = nstart + (rel if nsign > 0 else nlen - 1 - rel)

# --- read the original map and rebin into the new layout
hic = hicstraw.HiCFile(HIC)
mzd = hic.getMatrixZoomData("assembly", "assembly", "observed", "NONE", "BP", READ_RES)
nb = int(np.ceil(total_new / PLOT_RES))
N = np.zeros(nb * nb)
CHUNK = 20_000_000                         # read upper-triangle records in row chunks (a dense 315-Mb read crashes hicstraw)
for a in range(0, L0, CHUNK):
    b_ = min(a + CHUNK, L0)
    rec = mzd.getRecords(a, b_ - 1, a, L0 - 1)
    x = np.fromiter((r.binX for r in rec), dtype=np.int64, count=len(rec))
    y = np.fromiter((r.binY for r in rec), dtype=np.int64, count=len(rec))
    v = np.fromiter((r.counts for r in rec), dtype=np.float64, count=len(rec))
    del rec
    m = (x <= y) & (x >= a) & (x < b_)     # each upper-triangle pixel once
    x, y, v = x[m] // READ_RES, y[m] // READ_RES, v[m]
    px, py = new_pos[x], new_pos[y]
    ok = (px >= 0) & (py >= 0)
    bx, by, v = px[ok] // PLOT_RES, py[ok] // PLOT_RES, v[ok]
    N += np.bincount(bx * nb + by, weights=v, minlength=nb * nb)
    off = x[ok] != y[ok]                    # mirror off-diagonal pixels
    N += np.bincount(by[off] * nb + bx[off], weights=v[off], minlength=nb * nb)
N = N.reshape(nb, nb)
np.savez_compressed("contacts_chr1-39_250kb.npz", N=N, res=PLOT_RES,
                    names=[c for c, _, _ in chrom_bounds], starts=[s for _, s, _ in chrom_bounds], ends=[e for _, _, e in chrom_bounds])

# --- report: share of contacts within chromosomes, and cis contact density per chromosome
cid = np.zeros(nb, dtype=int)
for j, (cn, s, e) in enumerate(chrom_bounds):
    cid[s // PLOT_RES: int(np.ceil(e / PLOT_RES))] = j
same = cid[:, None] == cid[None, :]
print(f"contacts within the same chromosome: {100 * N[same].sum() / N.sum():.1f}%")
print("chromosome\tlength_Mb\tcis_contacts_per_Mb\ttrans_share_%")
for j, (cn, s, e) in enumerate(chrom_bounds):
    r = slice(s // PLOT_RES, int(np.ceil(e / PLOT_RES)))
    cis = N[r, r].sum() / 2; tot = N[r, :].sum()
    print(f"{cn}\t{(e - s) / 1e6:.2f}\t{cis / ((e - s) / 1e6):.0f}\t{100 * (1 - N[r, r].sum() / tot):.1f}")

# --- plot
fig, ax = plt.subplots(figsize=(7.2, 6.4))
vmax = np.quantile(N[N > 0], 0.98)
im = ax.imshow(np.log10(N + 1), cmap="Reds", vmin=0, vmax=np.log10(vmax + 1), interpolation="none",
               extent=(0, total_new / 1e6, total_new / 1e6, 0))
mids = [((s + e) / 2) / 1e6 for _, s, e in chrom_bounds]
labels = [cn.replace("HiC_scaffold_", "") for cn, _, _ in chrom_bounds]
for _, s, e in chrom_bounds[1:]:
    ax.axhline(s / 1e6, color="grey", lw=0.15); ax.axvline(s / 1e6, color="grey", lw=0.15)
ax.set_xticks(mids); ax.set_xticklabels(labels, fontsize=4.5, rotation=90)
ax.set_yticks(mids); ax.set_yticklabels(labels, fontsize=4.5)
ax.tick_params(length=1.5, width=0.4)
ax.set_xlabel("Chromosome (HiC_scaffold)", fontsize=8); ax.set_ylabel("Chromosome (HiC_scaffold)", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("log10(contacts + 1)", fontsize=7)
cb.ax.tick_params(labelsize=6)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"Figure1b_OmniC_final_assembly.{ext}", dpi=400)
print(f"saved Figure1b_OmniC_final_assembly.pdf/.png ({total_new / 1e6:.1f} Mb shown, {PLOT_RES // 1000} kb bins)")

# --- "clean" version: balanced matrix (iterative correction), log colour scale with everything at or
#     below the 99th percentile of between-chromosome contacts shown as white
good = N.sum(1) > 0.3 * np.median(N.sum(1)[N.sum(1) > 0])     # bins with too few reads are left white
B = N.copy(); B[~good] = 0; B[:, ~good] = 0
for _ in range(100):
    s = B.sum(1); s = s / np.mean(s[good]); s[~good] = 1; B = B / s[:, None] / s[None, :]
offdiag = ~np.eye(nb, dtype=bool)
lo = np.quantile(B[~same & (B > 0)], 0.99)
hi = np.quantile(B[same & offdiag & (B > 0)], 0.99)
fig, ax = plt.subplots(figsize=(7.2, 6.4))
cmap = matplotlib.colormaps["Reds"].copy(); cmap.set_bad("white")
im = ax.imshow(np.ma.masked_less_equal(np.log10(np.clip(B, lo, None)), np.log10(lo)), cmap=cmap,
               vmin=np.log10(lo), vmax=np.log10(hi), interpolation="none", extent=(0, total_new / 1e6, total_new / 1e6, 0))
ax.set_xticks(mids); ax.set_xticklabels(labels, fontsize=4.5, rotation=90)
ax.set_yticks(mids); ax.set_yticklabels(labels, fontsize=4.5)
ax.tick_params(length=1.5, width=0.4)
ax.set_xlabel("Chromosome (HiC_scaffold)", fontsize=8); ax.set_ylabel("Chromosome (HiC_scaffold)", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("log10(balanced contacts)", fontsize=7)
cb.ax.tick_params(labelsize=6)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"Figure1b_OmniC_final_assembly_clean.{ext}", dpi=400)
print("saved Figure1b_OmniC_final_assembly_clean.pdf/.png")

# --- middle version: observed (unbalanced) counts, log scale; full range within chromosomes,
#     between-chromosome contacts below their 75th percentile shown as white
TRANS_WHITE_Q = 0.75
lo = np.quantile(N[~same & (N > 0)], TRANS_WHITE_Q)
hi = np.quantile(N[same & offdiag & (N > 0)], 0.995)
fig, ax = plt.subplots(figsize=(7.2, 6.4))
im = ax.imshow(np.ma.masked_less_equal(np.log10(np.clip(N, lo, None)), np.log10(lo)), cmap=cmap,
               vmin=np.log10(lo), vmax=np.log10(hi), interpolation="none", extent=(0, total_new / 1e6, total_new / 1e6, 0))
ax.set_xticks(mids); ax.set_xticklabels(labels, fontsize=4.5, rotation=90)
ax.set_yticks(mids); ax.set_yticklabels(labels, fontsize=4.5)
ax.tick_params(length=1.5, width=0.4)
ax.set_xlabel("Chromosome (HiC_scaffold)", fontsize=8); ax.set_ylabel("Chromosome (HiC_scaffold)", fontsize=8)
cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("log10(contacts)", fontsize=7)
cb.ax.tick_params(labelsize=6)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"Figure1b_OmniC_final_assembly_observed.{ext}", dpi=400)
print(f"saved Figure1b_OmniC_final_assembly_observed.pdf/.png (white below {lo:.0f} contacts per bin pair)")
