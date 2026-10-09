"""Scan genomes for a satellite repeat: how much sequence sits in windows dense with its 12-mers.

Usage:  python satellite_scan.py satellite_296bp_consensus.fa genome1.fa[.gz] [genome2.fa.gz ...]

Every 12-mer of the monomer (circular, both strands) is looked up along each genome. In 2-kb windows
the share of positions that match is reported; random sequence gives ~0.004%, a homologous satellite
array gives several % even at ~75-80% identity (identical arrays reach 50-90%).
Needs numpy only.
"""
import sys, gzip
import numpy as np

K, WIN = 12, 2000
mono = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">")).upper()
rc = lambda s: s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
circ = mono + mono[:K - 1]
kms = {circ[i:i + K] for i in range(len(mono))} | {rc(circ)[i:i + K] for i in range(len(mono))}
sat = np.array(sorted(int("".join(str("ACGT".index(c)) for c in k), 4) for k in kms), dtype=np.uint64)
code = np.full(256, 255, np.uint8)
for i, b in enumerate("ACGT"):
    code[ord(b)] = i; code[ord(b.lower())] = i


def scan(path):
    op = gzip.open if path.endswith(".gz") else open
    dens, total = [], 0

    def do(seq):
        nonlocal total
        a = code[np.frombuffer(seq.encode(), np.uint8)]; n = len(a); total += n
        if n < WIN:
            return
        a2 = np.where(a == 255, 0, a).astype(np.uint64)
        h = np.zeros(n - K + 1, np.uint64)
        for j in range(K):
            h = h * np.uint64(4) + a2[j:n - K + 1 + j]
        m = np.isin(h, sat).astype(np.float32)
        nw = len(m) // WIN
        dens.append(m[:nw * WIN].reshape(nw, WIN).mean(1))

    buf = []
    for line in op(path, "rt"):
        if line.startswith(">"):
            if buf:
                do("".join(buf)); buf = []
        else:
            buf.append(line.strip())
    if buf:
        do("".join(buf))
    return (np.concatenate(dens) if dens else np.zeros(1)), total


print("genome\tsize_Mb\tMb_windows_>1%\tMb_windows_>10%\tmax_window_%")
for p in sys.argv[2:]:
    d, total = scan(p)
    print(f"{p}\t{total/1e6:.1f}\t{(d > 0.01).sum()*WIN/1e6:.2f}\t{(d > 0.10).sum()*WIN/1e6:.2f}\t{100*d.max():.1f}", flush=True)
