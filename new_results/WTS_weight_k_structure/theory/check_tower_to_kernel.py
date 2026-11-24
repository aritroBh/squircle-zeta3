"""LEMMA (tower -> kernel), exact check.  For every level k (2..14), every key (i,j) of degree d (1<=d<=k-1):
the kernel K^{(k)}_{ij} equals the real-pole part of
      -(i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2 - [j=i-1] i F(Gm_{ii} DL1),
where Gm_{ab} = (1-eval_b) o Int of the recorded +i principal part (orders >=2) H_{ab}^{(k)} of the parent work item
(and its complex conjugate at -i), i.e. the kernel is a k-INDEPENDENT Q(sqrt2,i)-linear image of the degree-(d+1) tower
coefficients; the +-i simple-pole part of that expression vanishes; the initial-integrand term eps^d R_ij only enters at k=d.
Also checks that the higher part (orders>=2 at +i) of the same expression is exactly the recorded tower coefficient of the key.
Uses tower_k.pkl (recorded by tower.py)."""
import sys, os, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from fastred import *

def q4(q): return (q.p0, q.p1, q.p2, q.p3)
def Q(t): return QSI(*t)

def Gm_from_tower(Hd):
    """Hd: {order r: 4-tuple} (+i coefficients, r>=2).  Build Gm = (1-eval_b) Int HigherPart, with -i part the conjugate."""
    L = []
    if Hd:
        top = max(Hd)
        L = [ZERO]*top
        for r, c in Hd.items(): L[r-1] = Q(c)
    H = pf_zero(); H[1] = list(L); H[2] = [e.conj_i() for e in L]
    G = pf_integrate_higher(H)
    C = pf_eval_at(G, BQ)
    return [G[0], G[1], G[2], -C]

def kernel_from_tower(T, k, key):
    i, j = key
    tot = pf_zero()
    if (i+1, j) in T:
        Gm = Gm_from_tower(T[(i+1, j)])
        tot = pf_add(tot, pf_scale(pf_mul_DL(Gm, DL1), Fr(-(i+1))))
    if i > j and (i, j+1) in T:
        Gm = Gm_from_tower(T[(i, j+1)])
        tot = pf_add(tot, pf_scale(pf_mul_DL(Gm, DL2), Fr(-(j+1))))
    if j == i-1 and (i, i) in T:
        Gm = Gm_from_tower(T[(i, i)])
        tot = pf_add(tot, pf_scale(pf_fold(pf_mul_DL(Gm, DL1)), Fr(-i)))
    return tot

allok = True
for fn in sorted(glob.glob(os.path.join(HERE, "tower_*.pkl")), key=lambda f: int(os.path.basename(f)[6:-4])):
    k = int(os.path.basename(fn)[6:-4])
    T = pickle.load(open(fn, 'rb'))
    ker = pickle.load(open(os.path.join(HERE, f"kern_{k}.pkl"), 'rb'))
    nchk = 0; bad = []
    for key in ker:
        i, j = key; d = i+j
        if d < 1 or d >= k: continue            # d=k: initial-integrand simple part contributes as well
        P = kernel_from_tower(T, k, key)
        simple, higher, const = pf_split(P)
        real = {(r.p0, r.p1): q4(c) for r, c in simple[0].items()}
        has_i = any(not e.is_zero() for e in simple[1]) or any(not e.is_zero() for e in simple[2]) or (not const.is_zero())
        # the higher part (orders >=2 at +i) of the same feed is the recorded tower coefficient of this key
        rec = {r: c for r, c in T.get(key, {}).items() if any(x != 0 for x in c)}
        got = {r: q4(e) for r, e in enumerate(higher[1], start=1) if not e.is_zero()}
        if got != rec:
            bad.append(("tower-mismatch", key))
        # canonical comparison: drop zero residues
        real = {r: c for r, c in real.items() if any(x != 0 for x in c)}
        want = {r: c for r, c in ker[key]['real'].items() if any(x != 0 for x in c)}
        ok = (real == want) and not has_i
        nchk += 1
        if not ok: bad.append(key)
    allok = allok and not bad
    print(f"k={k}: kernels of degree 1..{k-1} (keys checked: {nchk}) equal tower->kernel linear image, no +-i simple poles: {not bad}  bad={bad}", flush=True)
print("TOWER->KERNEL LEMMA verified for all levels above:", allok)
