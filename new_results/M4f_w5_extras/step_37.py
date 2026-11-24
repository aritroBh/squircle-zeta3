WORK='<project>/research/08_theory/M4f_w5_extras'
src = r'''"""M4f step 3b: exact atom-height optimizer over each subset's affine solution space.

For every closing subset S (subset_exact.pkl): v(t) = v0 + sum t_i n_i (146-dim Fraction,
k = 2|S| dims). Minimize atom height:
  score = (max|numerator|, lcm denominators, #atoms nonzero)
after display reductions Li5(b^2)->16(Li5(b)+Li5(-b)), Li5(-1)->-15/16 zeta(5).
Method: D-scaled integer LLL (gram exact) + Babai CVP + greedy integer/fractional descent.
Saves optimized.pkl: {S: (score, v)}.
"""
import json, pickle, sys, math, itertools, random
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
import nice_search as NS

mp.mp.dps = 100
def log(s): print(s, flush=True)

sol = NS.solve(list(NS.S0), 'opt')
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

def optimize(v0, NSimg, label='', exhaustive=False):
    k = len(NSimg)
    D = 1
    for a in v0: D = lcm(D, a.denominator)
    for nv in NSimg:
        for a in nv: D = lcm(D, a.denominator)
    w0 = [int(a*D) for a in v0]
    W = [[int(a*D) for a in nv] for nv in NSimg]
    if k == 0:
        return score(v0), v0, []
    red = flint.fmpz_mat(W).lll(gram='exact')
    B = [[int(red[i,j]) for j in range(NATOM)] for i in range(k)]
    # Babai nearest plane with mpmath GSO
    Bm = [[mp.mpf(x) for x in row] for row in B]
    # GSO
    n = k
    Bs = []; mu = [[mp.mpf(0)]*n for _ in range(n)]; norm = []
    for i in range(n):
        v = Bm[i][:]
        for j in range(i):
            mu[i][j] = mp.fsum(Bm[i][s]*Bs[j][s] for s in range(NATOM))/norm[j]
            for s in range(NATOM): v[s] -= mu[i][j]*Bs[j][s]
        Bs.append(v); norm.append(mp.fsum(x*x for x in v))
    t = [mp.mpf(0)]*n
    w = [mp.mpf(x) for x in w0]
    tt = [0]*n
    for i in range(n-1, -1, -1):
        c = mp.fsum(w[s]*Bs[i][s] for s in range(NATOM))/norm[i]
        ci = int(mp.nint(c))
        tt[i] = ci
        for s in range(NATOM): w[s] -= ci*Bm[i][s]
    def build(tt):
        v = [F(w0[s], D) for s in range(NATOM)]
        for i,ti in enumerate(tt):
            if ti:
                for s in range(NATOM):
                    if B[i][s]: v[s] += F(ti)*F(B[i][s], D)
        return v
    best_tt = list(tt); best_v = build(tt); best_s = score(best_v)
    # greedy descent on reduced coords, integer steps then small fractions
    improved = True
    rnd = 0
    while improved and rnd < 40:
        improved = False; rnd += 1
        for i in range(n):
            for step in [F(1), F(-1), F(2), F(-2), F(1,2), F(-1,2), F(1,3), F(-1,3), F(2,3), F(-2,3)]:
                cand = list(best_tt)
                cand[i] = [F(x) for x in cand]
                cand[i] = cand  # noop
                t2 = [F(x) for x in best_tt] ; t2[i] += step
                v2 = build_frac(t2, w0, B, D)
                s2 = score(v2)
                if s2 < best_s:
                    best_s, best_tt, best_v = s2, t2, v2
                    improved = True
    # pairwise random integer combos
    rng = random.Random(12345)
    for _ in range(400 if exhaustive else 120):
        t2 = [F(x) for x in best_tt]
        for _ in range(3):
            i = rng.randrange(n); t2[i] += rng.choice([1,-1,2,-2])
        v2 = build_frac(t2, w0, B, D); s2 = score(v2)
        if s2 < best_s:
            best_s, best_tt, best_v = s2, t2, v2
    return best_s, best_v, best_tt

def build_frac(tt, w0, B, D):
    v = [F(w0[s], D) for s in range(NATOM)]
    for i,ti in enumerate(tt):
        if ti:
            for s in range(NATOM):
                if B[i][s]: v[s] += ti*F(B[i][s], D)
    return v

def main():
    results = pickle.loads((HERE/'subset_exact.pkl').read_bytes())
    out = {}
    # optimize all; full set exhaustively
    for r in results:
        S = r['S']
        ex = len(S) == 7
        sc, bv, tt = optimize(r['v0'], r['NSimg'], label=str(S), exhaustive=ex)
        out[S] = (sc, bv)
        log(f'S={S} k={len(r["NSimg"])}: max|num|={sc[0]} lcm={sc[1]} terms={sc[2]}')
    pickle.dump(out, open(HERE/'optimized.pkl','wb'))
    best = sorted(out.items(), key=lambda kv: kv[1][0])[:10]
    log('TOP 10:')
    for S,(sc,_) in best:
        log(f'  S={S}: max|num|={sc[0]} lcm={sc[1]} terms={sc[2]}')

if __name__=='__main__':
    main()
'''
open(f'{WORK}/optimize.py','w').write(src)
print('written')
