"""w6 identification step 11: LLL relation lattice among the 12 integrals N and the 13 level-2 basis elements.
Two precisions: (250-run N, D=200) and (400-run N, D=300... using the 250 digits of S6 only up to 260 -> D=240).  Printed from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
from w6lll import lll_relations
for Nrun, D in ((250, 200), (400, 240)):
    mp.mp.dps = D + 20
    C = consts(); B = basis13(C)
    N = json.load(open(os.path.join(W, f"Nints_{Nrun}.json")))["N"]
    mp.mp.dps = D + 20
    names = [f"N[{k}]" for k in N] + [nm for nm, _ in B]
    vals = [mp.mpf(v) for v in N.values()] + [v for _, v in B]
    res = lll_relations(vals, D)
    print(f"--- N run {Nrun}, D={D}, {len(vals)} numbers, mp.mp.dps={mp.mp.dps}")
    k = 0
    for rel, r, last in res:
        if abs(last) > 10 ** 6 or abs(r) > mp.mpf(10) ** (-(D - 30)):
            break
        k += 1
        nz = [(names[i], a) for i, a in enumerate(rel) if a != 0]
        print(f" relation {k}: max|coeff|={max(abs(a) for a in rel)}  residual {mp.nstr(r,3)}")
        print("    ", nz)
    print(f"  number of relations found (residual < 1e-{D-30}): {k}")
