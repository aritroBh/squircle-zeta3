"""M4f step 3c: cascaded exact atom-height optimizer.

Affine space per subset S: v(t) = v0 + sum t_i n_i. Spaces nest (S' subset S).
Optimize bottom-up, seeding each S with the best point of its best optimized subset.
Descent: per-coordinate sweep over step sizes; then pairwise; score as in optimize.py.
Saves optimized2.pkl: {S: (score, v, t)}.
"""
import json, pickle, sys, math, random
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
import nice_search as NS

mp.mp.dps = 120
def log(s): print(s, flush=True)
sol = NS.solve(list(NS.S0), 'opt2')
names = sol['names']
ib2 = names.index('Li5((3-2*sqrt2))'); ibp = names.index('Li5((-1+1*sqrt2))')
ibm = names.index('Li5((1-1*sqrt2))'); im1 = names.index('Li5(-1)'); iz5 = names.index('z5')
NATOM = 146
PAIRS = [(i, i+65) for i in range(65)] + [(130,131)] + [(132+i,139+i) for i in range(7)]
def reduced(v):
    v = list(v)
    for off in (0,65):
        a = v[ib2+off]
        if a: v[ibp+off] += 16*a; v[ibm+off] += 16*a; v[ib2+off] = F(0)
        a = v[im1+off]
        if a: v[iz5+off] += F(-15,16)*a; v[im1+off] = F(0)
    return v
def score(v):
    v = reduced(v)
    nz = [a for a in v if a]
    terms = sum(1 for p,q in PAIRS if v[p] or v[q])
    return (max(abs(a.numerator) for a in nz), lcm(*[a.denominator for a in nz]), terms)

def build(tt, w0, B, D):
    v = [F(w0[s], D) for s in range(NATOM)]
    for i,ti in enumerate(tt):
        if ti:
            for s in range(NATOM):
                if B[i][s]: v[s] += ti*F(B[i][s], D)
    return v

class Space:
    def __init__(self, v0, NSimg):
        self.k = len(NSimg)
        D = 1
        for a in v0: D = lcm(D, a.denominator)
        for nv in NSimg:
            for a in nv: D = lcm(D, a.denominator)
        self.D = D
        self.w0 = [int(a*D) for a in v0]
        W = [[int(a*D) for a in nv] for nv in NSimg]
        if self.k:
            red = flint.fmpz_mat(W).lll(gram='exact')
            self.B = [[int(red[i,j]) for j in range(NATOM)] for i in range(self.k)]
        else:
            self.B = []
        self.W = W  # original (unreduced) for coordinate transforms
    def t_to_v(self, tt):
        return build(tt, self.w0, self.B, self.D)
    def coords_of(self, v):
        """Given v in the affine space, find tt in the REDUCED basis coords."""
        diff = [F(a) - F(w0s, self.D) for a, w0s in zip(v, self.w0)]
        # solve sum tt_i B_i = diff over Q (B rows independent)
        M = [[str(F(self.B[j][s], self.D)) for j in range(self.k)] + [str(diff[s])] for s in range(NATOM)]
        mat, rank = flint.fmpq_mat(M).rref()
        pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
        assert all(p is None or p < self.k for p in pivots), 'v not in space'
        tt = [F(0)]*self.k
        for i,p in enumerate(pivots):
            if p is not None and p < self.k: tt[p] = F(str(mat[i,self.k]))
        return tt
    def descend(self, tt, passes=30, steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2),F(1,3),F(-1,3),F(2,3),F(-2,3),F(1,6),F(-1,6))):
        best_v = self.t_to_v(tt); best_s = score(best_v); tt = list(tt)
        for _ in range(passes):
            improved = False
            for i in range(self.k):
                for st in steps:
                    t2 = list(tt); t2[i] += st
                    v2 = self.t_to_v(t2); s2 = score(v2)
                    if s2 < best_s: best_s, tt, best_v = s2, t2, v2; improved = True
            # pairwise
            for i in range(self.k):
                for j in range(i+1, self.k):
                    for si in (F(1),F(-1)):
                        for sj in (F(1),F(-1)):
                            t2 = list(tt); t2[i] += si; t2[j] += sj
                            v2 = self.t_to_v(t2); s2 = score(v2)
                            if s2 < best_s: best_s, tt, best_v = s2, t2, v2; improved = True
            if not improved: break
        return best_s, best_v, tt

def main():
    results = pickle.loads((HERE/'subset_exact.pkl').read_bytes())
    spaces = {r['S']: Space(r['v0'], r['NSimg']) for r in results}
    out = {}
    for S in sorted(spaces, key=lambda s:(len(s), s)):
        sp = spaces[S]
        # seed candidates: Babai-free zero seed, plus best point of best strict subset
        seeds = [[F(0)]*sp.k]
        subs = [s2 for s2 in out if set(s2) < set(S)]
        if subs:
            bestsub = min(subs, key=lambda s2: out[s2][0])
            seeds.append(sp.coords_of(out[bestsub][1]))
        best = None
        for sd in seeds:
            sc, bv, tt = sp.descend(sd)
            if best is None or sc < best[0]: best = (sc, bv, tt)
        out[S] = best
        log(f'S={S} k={sp.k}: max|num|={best[0][0]} lcm={best[0][1]} terms={best[0][2]}')
    pickle.dump(out, open(HERE/'optimized2.pkl','wb'))
    log('TOP 12:')
    for S,(sc,_,_) in sorted(out.items(), key=lambda kv: kv[1][0])[:12]:
        log(f'  S={S}: max|num|={sc[0]} lcm={sc[1]} terms={sc[2]}')

if __name__=='__main__':
    main()
