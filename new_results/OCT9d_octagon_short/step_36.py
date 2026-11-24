import sys, time, itertools, pickle
sys.path.insert(0, ".")
from oct9d_common import *
t0=time.time()
TW = TW3()
def colbasis(M):
    rref, pivs = M.rref(); return [M[:,i] for i in pivs]

# Li3 symbol columns in V (27-dim) for all Gargs
Gsym = {g: col_Li3(sp.Rational(g), SS235) for g in Gargs}
# product columns, nice order: w small height, p in (2,3,5)
EP = set()
for e in (1,2,3):
    for z in edge_points(e): EP.add(z)
ORB = set()
for z in list(EP): ORB |= orbit(z)
ORBSM = set()
for z in list(SMALL235): ORBSM |= orbit(z)
def nice(g):
    g = Fr(g)
    if g in EP: return 0
    if g in ORB: return 1
    if g in SMALL235: return 2
    if g in ORBSM: return 3
    return 4
worder = sorted(Gargs, key=lambda g: (nice(g), max(abs(Fr(g).numerator), abs(Fr(g).denominator)), str(g)))
PB = []  # (w,p, col)
cur = sp.zeros(27,0)
for w in worder:
    for p in (2,3,5):
        c = col_prod(sp.Rational(w), p, SS235)
        if sp.Matrix.hstack(cur, c).rank() > cur.shape[1]:
            PB.append((w, p, c)); cur = sp.Matrix.hstack(cur, c)
        if cur.shape[1] == 19: break
    if cur.shape[1] == 19: break
print("product basis (19):", [(str(w), p) for w,p,c in PB], f"({time.time()-t0:.1f}s)")
Pmat = sp.Matrix.hstack(*[c for w,p,c in PB])

cand_bases = [
    ['1/3', '2/3', '2/5', '-1/5', '1/10', '-5/3', '5', '3/5'],
    ['10/9', '-2/3', '-1/2', '1/3', '-5/3', '-1/5', '5/4', '2/5'],
    ['1/3', '2/5', '5/8', '-5/3', '2/3', '9/5', '-2/3', '4/5'],
]
for sel in cand_bases:
    selq = [sp.Rational(g) for g in sel]
    M = sp.Matrix.hstack(*[Gsym[g] for g in selq], Pmat)
    print("basis", sel, "det:", M.det())
    for e in (1,2,3):
        x = M.LUsolve(sp.Matrix(TW[e]))
        cl = x[:8]; cp = x[8:]
        print(f"  e{e} Li3 coeffs:", [str(c) for c in cl])
        print(f"  e{e} prod coeffs:", [str(c) for c in cp])
