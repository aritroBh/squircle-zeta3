"""verify step: exact-rational check (own parser, own finite differences) that the residues of the
weight-(d+1) kernels K_ij (i+j=d) of the fast-reducer output finals_k.json are polynomials in k.
For d = 5 (weight 6): k = 6..15 must be degree <= 4 (5th exact differences zero), and k = 5 must break it.
Also scans d = 2,3,4,6,7 for the minimal exact degree on k >= d+1.
usage: verify_exact_poly.py
"""
import sys, os, json
sys.dont_write_bytecode = True
W = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, W)
from fractions import Fraction as Fr
import verify_lib as VL

KS = list(range(3, 16))
F = {k: VL.load_finals(os.path.join(W, f"finals_{k}.json")) for k in KS}

def diffs(seq, n):
    cur = list(seq)
    for _ in range(n):
        cur = [cur[i + 1] - cur[i] for i in range(len(cur) - 1)]
    return cur

def comps(d):
    """(key,pole,component) triples occurring for degree d at any k"""
    out = set()
    for k in KS:
        if k < d:
            continue
        for key, poles in F[k].items():
            if sum(key) != d:
                continue
            for lab in poles:
                out.add((key, lab))
    return sorted(out)

def val(k, key, lab, comp):
    p = F[k].get(key, {}).get(lab)
    return Fr(0) if p is None else p[comp]

for d in (1, 2, 3, 4, 5, 6, 7):
    kmin = d + 1
    ks = [k for k in range(max(kmin, 3), 16)]
    allc = comps(d)
    mind = 0
    detail = []
    for key, lab in allc:
        for comp in (0, 1):
            ys = [val(k, key, lab, comp) for k in ks]
            if all(y == 0 for y in ys):
                continue
            # minimal exact degree on ks
            dg = None
            for n in range(1, len(ys)):
                if all(x == 0 for x in diffs(ys, n)):
                    dg = n - 1
                    break
            detail.append((key, lab, comp, dg))
            mind = max(mind, dg if dg is not None else 99)
    n_nonzero = len(detail)
    print(f"d={d} (weight {d+1}): k={ks[0]}..{ks[-1]} ({len(ks)} pts); nonzero (key,pole,component) = {n_nonzero}; "
          f"max exact degree = {mind}; polynomial for every component: {all(x[3] is not None for x in detail)}")
    if d == 5:
        for (key, lab, comp, dg) in detail:
            ys6 = [val(k, key, lab, comp) for k in range(6, 16)]
            ys5 = [val(k, key, lab, comp) for k in range(5, 16)]
            # does including k=5 keep the degree <= dg?
            ok5 = all(x == 0 for x in diffs(ys5, dg + 1)) if len(ys5) > dg + 1 else None
            print(f"    key {key} pole {lab:>2} comp {'p0' if comp==0 else 'p1'}: exact degree on k=6..15 = {dg}; "
                  f"still degree {dg} when k=5 included: {ok5}")
        print("    poles present at k>=6 for d=5:", sorted({lab for (key, lab, comp, dg) in detail}), "(b absent => all keys converge)")
        # k = 5 itself: which keys/poles
        print("    k=5, d=5 keys:", {key: sorted(p) for key, p in F[5].items() if sum(key) == 5})
        print("    k=6, d=5 keys:", {key: sorted(p) for key, p in F[6].items() if sum(key) == 5})
