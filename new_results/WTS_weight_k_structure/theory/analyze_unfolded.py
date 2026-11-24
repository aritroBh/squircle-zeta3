"""Side experiment: do the kernels of the UNFOLDED reduction (no v->1/v fold; keys (i,j) with all i,j) obey the same
polynomial-in-k law?  Exact finite differences on ukern_k.pkl for k>=d+1.  Reports minimal degree per key vs d-1."""
import sys, os, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
from fractions import Fraction as Fr
U = {int(os.path.basename(f)[6:-4]): pickle.load(open(f, 'rb')) for f in glob.glob(os.path.join(HERE, "ukern_*.pkl"))}
ks = sorted(U); K = ks[-1]; print("unfolded levels:", ks)
def fd(seq, m):
    s = list(seq)
    for _ in range(m): s = [s[i+1]-s[i] for i in range(len(s)-1)]
    return s
keys = sorted({key for k in ks for key in U[k] if sum(key) >= 1}, key=lambda kk: (sum(kk), kk))
for key in keys:
    d = sum(key)
    kr = [k for k in ks if k >= d+1]
    if len(kr) < d+2: continue
    poles = sorted({r for k in kr for r in U[k].get(key, {'real': {}})['real']})
    best = -1
    for r in poles:
        for comp in (0, 1):
            seq = [U[k].get(key, {'real': {}})['real'].get(r, (Fr(0),)*4)[comp] for k in kr]
            mdeg = None
            for m in range(len(seq)):
                if all(x == 0 for x in fd(seq, m)): mdeg = m-1; break
            best = max(best, 99 if mdeg is None else mdeg)
    print(f"key {key} d={d}: minimal poly degree in k (k>=d+1) = {best if best < 99 else 'NOT polynomial in tested range'};  conjectured bound d-1 = {d-1};  points {len(kr)}")
