"""Re-run the fast reducer for level k and RECORD, for every key (i,j) (i>=j) the full work item W_ij^{(k)}:
   its higher part at +i (coefficients e_r, r>=2, in Q(sqrt2,i)) and the simple part's real residues.
usage: python tower.py K1 K2 -> tower_k.pkl for K1<=k<=K2 ; each: {'H':{key:{r:(4 Fractions)}}, 'SIMPLE':{key:{pole:(c0..c3)}} }
Same arithmetic as fastred.reduce_k_fast (code copied line by line; only recording added)."""
import sys, os, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from fastred import *
from buildintegrand import build_integrand_pf

def q4(q): return (q.p0, q.p1, q.p2, q.p3)

def reduce_k_tower(k, integrand_pf):
    rec_H = {}; rec_W = {}
    work = {}
    for key, R in integrand_pf.items():
        kk, RR = fold_key(key, R)
        work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    final = {}
    for deg in range(k, -1, -1):
        for key in sorted([kk for kk in work if sum(kk) == deg], reverse=True):
            R = work.pop(key)
            if deg == 0:
                final[key] = R
                continue
            simple, higher, const = pf_split(R)
            assert const.is_zero()
            rec_H[key] = {r: q4(e) for r, e in enumerate(higher[1], start=1) if not e.is_zero()}
            if simple[0] or simple[1] or simple[2]:
                final[key] = pf_add(final.get(key, pf_zero()), simple)
            if not higher[1] and not higher[2]:
                continue
            G = pf_integrate_higher(higher)
            C = pf_eval_at(G, BQ)
            Gm = [G[0], G[1], G[2], -C]
            i, j = key
            if i >= 1:
                kk, RR = fold_key((i-1, j), pf_scale(pf_mul_DL(Gm, DL1), Fr(-i)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
            if j >= 1:
                kk, RR = fold_key((i, j-1), pf_scale(pf_mul_DL(Gm, DL2), Fr(-j)))
                work[kk] = pf_add(work.get(kk, pf_zero()), RR)
    return final, rec_H

if __name__ == "__main__":
    k1, k2 = int(sys.argv[1]), int(sys.argv[2])
    for k in range(k1, k2+1):
        ipf = pickle.load(open(os.path.join(HERE, f"binteg_{k}.pkl"), 'rb'))
        t0 = time.time()
        final, H = reduce_k_tower(k, ipf)
        # cross-check against the saved kernels (same code path => identical)
        saved = pickle.load(open(os.path.join(HERE, f"kern_{k}.pkl"), 'rb'))
        for key, P in final.items():
            sp_ = saved[key]
            assert {(r.p0, r.p1): q4(c) for r, c in P[0].items()} == sp_['real'], ("kernel mismatch vs saved", k, key)
        assert set(final) == set(saved)
        pickle.dump(H, open(os.path.join(HERE, f"tower_{k}.pkl"), 'wb'))
        print(f"k={k}: tower recorded ({len(H)} keys) in {time.time()-t0:.1f}s; kernels identical to saved kern_{k}.pkl", flush=True)
