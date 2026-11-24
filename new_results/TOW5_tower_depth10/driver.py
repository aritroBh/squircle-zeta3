
"""harvest + tower-law test driver (Task B)."""
import os, sys, pickle, time
from fractions import Fraction as Fr
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import trunc_engine as te
from fastred import QSI, ZERO

def q4(q): return (q.p0, q.p1, q.p2, q.p3)
def unq4(t): return QSI(*t)

def harvest(D, K, outpath, k_start=2):
    if os.path.exists(outpath):
        H = pickle.load(open(outpath,'rb'))
        done = H.pop('_meta', {}).get('k_done', 1)
        k_start = max(k_start, done+1)
        H = {kk: vv for kk, vv in H.items() if kk != '_meta'}
    else:
        H = {}
    S2 = te.stirling_table(K+1)
    for k in range(k_start, K+1):
        t0 = time.time()
        rho_k = te.rho_level(k, D, S2)
        tower_k = te.build_tower_level(k, D, rho_k)
        for n in range(2, k+1):
            for key in te.keys_of_degree(n):
                vals = tower_k[key]
                if key not in H: H[key] = {}
                for delta in range(0, min(D, 2*n-4)+1):
                    H[key].setdefault(delta, {})[k] = q4(vals[delta])
        out = dict(H); out['_meta'] = {'k_done': k, 'D': D, 'K': K}
        pickle.dump(out, open(outpath,'wb'))
        print(f"  level k={k} harvested ({time.time()-t0:.1f}s)", flush=True)
    return H

def fd_test(ys, P):
    m = len(ys)-1
    cur = list(ys)
    cols = [cur]
    for j in range(1, m+1):
        cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
        cols.append(cur)
    if m >= P+1:
        checks = cols[P+1]
        is_poly = all(c.is_zero() for c in checks)
        nchecks = len(checks)
    else:
        is_poly = None; nchecks = 0
    deg_exact = (m >= P) and (not cols[P][0].is_zero())
    return is_poly, deg_exact, m+1, nchecks

def test_tower_law(H, D, N, K):
    cells = []
    P_of = lambda delta: delta//2
    for key in sorted(H.keys(), key=lambda k:(k[0]+k[1], -k[0])):
        p, q = key; n = p+q
        if n > N: continue
        diag = (p == q)
        k0 = n+1 if diag else n
        for delta in range(0, min(D, 2*n-4)+1):
            kd = H[key].get(delta, {})
            ks = sorted(k for k in kd.keys() if k >= k0)
            ys = [unq4(kd[k]) for k in ks]
            P = P_of(delta)
            is_poly, deg_exact, npts, nchecks = fd_test(ys, P)
            genuine = (nchecks >= 1)
            cells.append(dict(key=key, n=n, delta=delta, P=P, diag=diag, k0=k0,
                              kmax=ks[-1] if ks else None, npts=npts, nchecks=nchecks,
                              is_poly=is_poly, deg_exact=deg_exact, genuine=genuine))
    return cells
