
import os, pickle
from fractions import Fraction as Fr
BASE = '<project>'
WT2B = os.path.join(BASE, 'research/08_theory/WT2b_stabilization')
rec = pickle.load(open(os.path.join(WT2B,'record_W2b.pkl'),'rb'))

def fit_poly_fr(vals_ks, deg):
    """vals_ks: dict k->Fr. Fit polynomial in k of given degree via finite differences / lagrange; return coeffs in binomial basis."""
    ks = sorted(vals_ks)
    ys = [vals_ks[k] for k in ks]
    # finite differences at k0=ks[0]
    fd = []
    cur = ys[:]
    for r in range(deg+2):
        fd.append(cur[0])
        cur = [cur[i+1]-cur[i] for i in range(len(cur)-1)]
    return ks[0], fd  # value expressed as sum_r fd[r] * binom(k-k0, r)

for key in [(4,0),(3,1),(5,0),(4,1),(3,2)]:
    d = sum(key)
    print(f"===== Gm {key} (d={d}): binomial-basis finite differences (k0, [c0,c1,c2,...]) =====")
    ks = sorted(k for (k,kk) in rec if kk==key)
    for m in [1,2,3,4,5]:
        vals = {}
        for k in ks:
            G = rec[(k,key)]['Gm']
            if m <= len(G[1]):
                vals[k] = G[1][m-1]
        if not vals: continue
        for ci, nm in [(0,'p0'),(1,'p1'),(2,'p2'),(3,'p3')]:
            vks = {k: getattr(c, nm) for k, c in vals.items()}
            if all(x==0 for x in vks.values()): continue
            k0, fd = fit_poly_fr(vks, 5)
            print(f"  ord{m} {nm}: k0={k0} fd={fd}")
    # const
    vals = {k: rec[(k,key)]['Gm'][3] for k in ks}
    for ci, nm in [(0,'p0'),(1,'p1'),(2,'p2'),(3,'p3')]:
        vks = {k: getattr(c, nm) for k, c in vals.items()}
        if all(x==0 for x in vks.values()): continue
        k0, fd = fit_poly_fr(vks, 5)
        print(f"  const {nm}: k0={k0} fd={fd}")
