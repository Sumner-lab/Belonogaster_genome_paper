"""Figure 1: (a) B. juncea photo, (b) Omni-C contact map of the 39 chromosome-scale scaffolds.

Both panels are drawn as squares of the same size. Panel b uses the same colour scaling as
HiC_replot/replot_omnic_final.py ("observed" version): log10 observed contacts in 250-kb bins,
between-chromosome bin pairs below the 75th percentile shown white.

Usage: ~/venvs/rxiv/bin/python3 make_figure1.py   (needs numpy, matplotlib, pillow)
Writes Figure1.pdf/.png/.tif (axis labels 4.5 pt, all 39 labelled) and Figure1_largefont.* (axis labels 6.5 pt).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

PHOTO = "photo_square.tif"  # DSC_0743.NEF, orientation tag stripped (pixels are stored heads-up), square crop, brightened (see README)
NPZ = "contacts_chr1-39_250kb.npz"
TRANS_WHITE_Q = 0.75

d = np.load(NPZ)
N, res = d["N"], int(d["res"])
bounds = list(zip([str(n) for n in d["names"]], d["starts"].astype(int), d["ends"].astype(int)))
total = bounds[-1][2]
nb = N.shape[0]
cid = np.zeros(nb, dtype=int)
for j, (_, s, e) in enumerate(bounds):
    cid[s // res: int(np.ceil(e / res))] = j
same = cid[:, None] == cid[None, :]
offdiag = ~np.eye(nb, dtype=bool)
lo = np.quantile(N[~same & (N > 0)], TRANS_WHITE_Q)
hi = np.quantile(N[same & offdiag & (N > 0)], 0.995)
cmap = matplotlib.colormaps["Reds"].copy()
cmap.set_bad("white")
M = np.ma.masked_less_equal(np.log10(np.clip(N, lo, None)), np.log10(lo))

mids = np.array([(s + e) / 2 / 1e6 for _, s, e in bounds])
labels = [n.replace("HiC_scaffold_", "") for n, _, _ in bounds]
photo = Image.open(PHOTO)

MM = 1 / 25.4
FIG_W, FIG_H = 180 * MM, 93 * MM     # full page width
SIDE = 72 * MM                       # both panels are SIDE x SIDE


def thin_labels(pos_mm, min_gap_mm):
    """Keep a tick label only if it is at least min_gap_mm from the previous kept label."""
    keep, last = [], -1e9
    for i, p in enumerate(pos_mm):
        if p - last >= min_gap_mm:
            keep.append(i)
            last = p
    return keep


def draw(tick_pt, out, label_all=False):
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    bottom = 12 * MM / FIG_H
    w, h = SIDE / FIG_W, SIDE / FIG_H
    ax_a = fig.add_axes([3 * MM / FIG_W, bottom, w, h])
    ax_b = fig.add_axes([(3 + 72 + 16) * MM / FIG_W, bottom, w, h])
    cax = fig.add_axes([(3 + 72 + 16 + 72 + 2.5) * MM / FIG_W, bottom + h * 0.15, 2.5 * MM / FIG_W, h * 0.7])

    ax_a.imshow(photo)
    ax_a.set_axis_off()

    im = ax_b.imshow(M, cmap=cmap, vmin=np.log10(lo), vmax=np.log10(hi), interpolation="none",
                     extent=(0, total / 1e6, total / 1e6, 0))
    # label every chromosome where there is room; skip some of the small ones at the end otherwise
    pos_mm = mids / (total / 1e6) * SIDE / MM
    keep = list(range(len(labels))) if label_all else thin_labels(pos_mm, tick_pt * 0.3528 * 1.05)
    ax_b.set_xticks(mids[keep]); ax_b.set_xticklabels([labels[i] for i in keep], fontsize=tick_pt, rotation=90)
    ax_b.set_yticks(mids[keep]); ax_b.set_yticklabels([labels[i] for i in keep], fontsize=tick_pt)
    ax_b.tick_params(length=1.5, width=0.4, pad=1)
    if label_all:
        # where neighbouring labels would overlap, stagger them: every other one gets a longer tick
        min_gap = tick_pt * 0.3528 * 1.0   # one em between label centres
        flip = False
        for k, i in enumerate(keep):
            crowded = (k > 0 and pos_mm[i] - pos_mm[keep[k - 1]] < min_gap) or \
                      (k + 1 < len(keep) and pos_mm[keep[k + 1]] - pos_mm[i] < min_gap)
            flip = (not flip) if crowded else False
            if crowded and flip:
                for t in (ax_b.xaxis.get_major_ticks()[k], ax_b.yaxis.get_major_ticks()[k]):
                    t.tick1line.set_markersize(1.5 + tick_pt * 1.3)
                    t.set_pad(1 + tick_pt * 1.3)
    for sp in ax_b.spines.values():
        sp.set_linewidth(0.5)
    ax_b.set_xlabel("Chromosome", fontsize=tick_pt + 2.5, labelpad=2)
    ax_b.set_ylabel("Chromosome", fontsize=tick_pt + 2.5, labelpad=2)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("log10(contacts)", fontsize=tick_pt + 1)
    cb.ax.tick_params(labelsize=tick_pt, length=1.5, width=0.4)
    cb.outline.set_linewidth(0.4)

    for lab, ax in (("a", ax_a), ("b", ax_b)):
        x0 = ax.get_position().x0
        dx = 0 if lab == "a" else -9 * MM / FIG_W
        fig.text(x0 + dx, bottom + h + 1.5 * MM / FIG_H, lab, fontsize=12, fontweight="bold", va="bottom")

    for ext in ("pdf", "png", "tif"):
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tif" else {}
        fig.savefig(f"{out}.{ext}", dpi=600, **kw)
    plt.close(fig)
    print(f"{out}: {len(keep)} of {len(labels)} chromosome labels at {tick_pt} pt")


draw(4.5, "Figure1", label_all=True)
draw(6.5, "Figure1_largefont")
