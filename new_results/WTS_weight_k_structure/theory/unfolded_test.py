"""Does the v->1/v fold matter?  Reduce the UNFOLDED integrand (all keys (i,j), i+j=k, no fold, no symmetrisation)
with the same exact algorithm (principal parts at +-i integrated, normalised at b, DL1/DL2 feeds), then FOLD the resulting
kernels:  Kfold_(i,j) = K_(i,j) + F(K_(j,i))  (i>j),  Kfold_(i,i) = K_(i,i),  and compare with the fast folded reducer
(kern_k.pkl).   Exact (Q(sqrt2,i)) equality test per key.  usage: python unfolded_test.py K1 K2"""
import sys, os, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from fastred import *
from genintegrand import S_poly
from buildintegrand import sqrt_part_PF, AB_part_PF

def q4(q): return (q.p0, q.p1, q.p2, q.p3)
def plain(P):
    return {'real': {(r.p0, r.p1): q4(c) for r, c in P[0].items()},
            'ipos': [q4(e) for e in P[1]], 'ineg': [q4(e) for e in P[2]], 'const': q4(P[3])}

def unfolded_integrand(k):
    out = {}
    for m in range(0, k+1):
        out[(m, k-m)] = sqrt_part_PF(k, S_poly(k, m))
    A = AB_part_PF(k)
    out[(k, 0)] = pf_add(out[(k, 0)], A)            # A-part (l1^k)
    out[(0, k)] = pf_add(out[(0, k)], pf_fold(A))   # B-part (l2^k) = fold of the A-form
    return out

def reduce_unfolded(k, ipf):
    work = {key: R for key, R in ipf.items()}
    final = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
            R = work.pop(key)
            if deg == 0:
                final[key] = R; continue
            simple, higher, const = pf_split(R)
            assert const.is_zero(), ("polynomial part", key)
            if simple[0] or simple[1] or simple[2]:
                final[key] = pf_add(final.get(key, pf_zero()), simple)
            if not higher[1] and not higher[2]:
                continue
            G = pf_integrate_higher(higher)
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            i, j = key
            if i >= 1:
                kk = (i-1, j); RR = pf_scale(pf_mul_DL(Gm, DL1), Fr(-i))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk = (i, j-1); RR = pf_scale(pf_mul_DL(Gm, DL2), Fr(-j))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    return final

def pf_is_zero(P):
    return (not P[0]) and all(e.is_zero() for e in P[1]) and all(e.is_zero() for e in P[2]) and P[3].is_zero()

def from_plain(d):
    P = pf_zero()
    for (r0, r1), c in d['real'].items():
        P[0][QSI(r0, r1)] = QSI(*c)
    P[1] = [QSI(*e) for e in d['ipos']]; P[2] = [QSI(*e) for e in d['ineg']]; P[3] = QSI(*d['const'])
    return P

k1, k2 = int(sys.argv[1]), int(sys.argv[2])
for k in range(k1, k2+1):
    t0 = time.time()
    ipf = unfolded_integrand(k)
    fin = reduce_unfolded(k, ipf)
    for key, P in fin.items():
        if sum(key) >= 1:
            assert not P[1] and not P[2], ("pole at +-i survives in unfolded reduction", k, key)
    pickle.dump({key: plain(P) for key, P in fin.items()}, open(os.path.join(HERE, f"ukern_{k}.pkl"), 'wb'))
    saved = pickle.load(open(os.path.join(HERE, f"kern_{k}.pkl"), 'rb'))
    allok = True; report = []
    for key in sorted(saved, key=lambda kk: (sum(kk), -kk[0])):
        i, j = key
        if sum(key) == 0: continue
        Ku = fin.get(key, pf_zero())
        if i > j:
            Kt = fin.get((j, i), pf_zero())
            Kf = pf_add(Ku, pf_fold(Kt)) if not pf_is_zero(Kt) else Ku
        else:
            Kf = Ku
        Ks = from_plain(saved[key])
        same = pf_is_zero(pf_sub(Kf, Ks))
        # also is the unfolded kernel by itself symmetric-compatible? record whether folded-sum matches
        allok = allok and same
        report.append((key, same))
    bad = [key for key, s in report if not s]
    print(f"k={k}: unfolded-then-folded kernels == folded reducer for all keys? {allok}; mismatching keys: {bad}; {time.time()-t0:.1f}s", flush=True)
