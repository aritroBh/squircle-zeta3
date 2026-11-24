WORK='<project>/research/08_theory/M4f_w5_extras'
src = r'''"""M4f step 3d: heavy optimization of the full 14-dim affine space (and (1..6)).
1. Iteratively reweighted LLL-CVP (target the inf-norm).
2. Slot-killing: exactly zero the current worst slot via the freedom, re-descend.
3. Random-restart descent.
Score = (max|num|, lcm, terms) after display reductions. Saves optimized3.pkl.
"""
import pickle, sys, random
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
import nice_search as NS

mp.mp.dps = 150
def log(s): print(s, flush=True)
sol = NS.solve(list(NS.S0), 'opt3')
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

class Space:
    def __init__(self, v0, NSimg):
        self.k = len(NSimg)
        D = 1
        for a in v0: D = lcm(D, a.denominator)
        for nv in NSimg:
            for a in nv: D = lcm(D, a.denominator)
        self.D = D
        self.w0 = [int(a*D) for a in v0]
        self.W = [[int(a*D) for a in nv] for nv in NSimg]  # original integer basis
    def lattice(self, weights=None):
        W = self.W
        if weights:
            W = [[int(round(x*w)) for x,w in zip(row,weights)] for row in W]
        red = flint.fmpz_mat(W).lll(gram='exact')
        return [[int(red[i,j]) for j in range(red.ncols())] for i in range(red.nrows())]
    def t_to_v(self, tt, B):
        v = [F(w0s, self.D) for w0s in self.w0]
        for i,ti in enumerate(tt):
            if ti:
                for s in range(NATOM):
                    if B[i][s]: v[s] += ti*F(B[i][s], self.D)
        return v

def babai(w0, B, NATOM):
    n = len(B)
    Bm = [[mp.mpf(x) for x in row] for row in B]
    Bs = []; mu = []; norm = []
    for i in range(n):
        v = Bm[i][:]
        mui = []
        for j in range(i):
            m = mp.fsum(Bm[i][s]*Bs[j][s] for s in range(len(v)))/norm[j]
            mui.append(m)
            for s in range(len(v)): v[s] -= m*Bs[j][s]
        Bs.append(v); mu.append(mui); norm.append(mp.fsum(x*x for x in v))
    w = [mp.mpf(x) for x in w0]
    tt = [0]*n
    for i in range(n-1, -1, -1):
        c = mp.fsum(w[s]*Bs[i][s] for s in range(NATOM))/norm[i]
        ci = int(mp.nint(c))
        tt[i] = ci
        for s in range(NATOM): w[s] -= ci*Bm[i][s]
    return tt

def descend(sp, B, tt, passes=40,
            steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2),F(1,3),F(-1,3),F(2,3),F(-2,3),F(1,6),F(-1,6))):
    best_v = sp.t_to_v(tt, B); best_s = score(best_v); tt = list(tt)
    for _ in range(passes):
        improved = False
        for i in range(sp.k):
            for st in steps:
                t2 = list(tt); t2[i] += st
                v2 = sp.t_to_v(t2, B); s2 = score(v2)
                if s2 < best_s: best_s, tt, best_v = s2, t2, v2; improved = True
        for i in range(sp.k):
            for j in range(i+1, sp.k):
                for si in (F(1),F(-1)):
                    for sj in (F(1),F(-1)):
                        t2 = list(tt); t2[i] += si; t2[j] += sj
                        v2 = sp.t_to_v(t2, B); s2 = score(v2)
                        if s2 < best_s: best_s, tt, best_v = s2, t2, v2; improved = True
        if not improved: break
    return best_s, best_v, tt

def main():
    results = pickle.loads((HERE/'subset_exact.pkl').read_bytes())
    out2 = pickle.loads((HERE/'optimized2.pkl').read_bytes())
    Sfull = (0,1,2,3,4,5,6)
    r = [r for r in results if r['S']==Sfull][0]
    sp = Space(r['v0'], r['NSimg'])
    B = sp.lattice()
    # start from optimized2 best point
    best_sc, best_v, _ = out2[Sfull]
    log(f'start: {best_sc}')
    # coords of best_v in B
    diff = [a - F(w0s, sp.D) for a, w0s in zip(best_v, sp.w0)]
    M = [[str(F(B[j][s], sp.D)) for j in range(sp.k)] + [str(diff[s])] for s in range(NATOM)]
    mat, rank = flint.fmpq_mat(M).rref()
    pivots = [next((j for j in range(mat.ncols()) if mat[j,j] if False else 0), None) for i in range(0)]
    pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
    assert all(p is None or p < sp.k for p in pivots)
    tt0 = [F(0)]*sp.k
    for i,p in enumerate(pivots):
        if p is not None and p < sp.k: tt0[p] = F(str(mat[i,sp.k]))
    # iteratively reweighted CVP rounds
    best_tt = tt0
    for rnd in range(6):
        v = sp.t_to_v(best_tt, B)
        sc = score(v)
        weights = [1.0 + (abs(float(a))* sp.D / sc[0])**4 for a in v]
        Ww = [[int(round(x*w)) for x,w in zip(row,weights)] for row in sp.W]
        red = flint.fmpz_mat(Ww).lll(gram='exact')
        Bw = [[int(red[i,j]) for j in range(red.ncols())] for i in range(red.nrows())]
        # unweight the basis: Bw rows are weighted combos; recover unweighted via weights
        Bunw = [[F(Bw[i][s])/F(int(round(weights[s]))) for s in range(NATOM)] for i in range(sp.k)]
        # re-derive integer unweighted basis by solving: Bunw ~ B combos; just use Babai in weighted metric:
        w0w = [int(round(x*w)) for x,w in zip(sp.w0, weights)]
        ttw = babai(w0w, Bw, NATOM)
        # convert weighted combo to actual vector: v = w0w - sum ttw_i Bw_i, then unweight
        v2 = [F(w0w[s]) for s in range(NATOM)]
        for i,ti in enumerate(ttw):
            for s in range(NATOM):
                if Bw[i][s]: v2[s] -= ti*Bw[i][s]
        v2 = [F(v2[s])/F(int(round(weights[s])))/sp.D for s in range(NATOM)]
        sc2, bv2, tt2 = descend_from_v(sp, B, v2)
        log(f'reweighted round {rnd}: {sc2}')
        if sc2 < score(sp.t_to_v(best_tt, B)):
            # recover tt in B coords
            best_tt = coords_of(sp, B, bv2)
            best_v = bv2
    sc, bv, tt = descend(sp, B, best_tt)
    log(f'after reweight+descend: {sc}')
    # random restarts
    rng = random.Random(777)
    for rs in range(120):
        t2 = [x + F(rng.choice([1,-1,2,-2]))*rng.choice([F(1),F(1,2),F(1,3)]) for x in tt]
        sc2, bv2, tt2 = descend(sp, B, t2, passes=6,
                                steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2)))
        if sc2 < sc:
            sc, bv, tt = sc2, bv2, tt2
            log(f'restart {rs}: improved {sc}')
    log(f'FINAL full-space: {sc}')
    pickle.dump(dict(S=Sfull, score=sc, v=bv), open(HERE/'optimized3.pkl','wb'))

def coords_of(sp, B, v):
    diff = [a - F(w0s, sp.D) for a, w0s in zip(v, sp.w0)]
    M = [[str(F(B[j][s], sp.D)) for j in range(sp.k)] + [str(diff[s])] for s in range(NATOM)]
    mat, rank = flint.fmpq_mat(M).rref()
    pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
    tt = [F(0)]*sp.k
    for i,p in enumerate(pivots):
        if p is not None and p < sp.k: tt[p] = F(str(mat[i,sp.k]))
    return tt

def descend_from_v(sp, B, v):
    tt = coords_of(sp, B, v)
    return descend(sp, B, tt)

if __name__=='__main__':
    main()
'''
open(f'{WORK}/optimize3.py','w').write(src)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','optimize3.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=2700)
print('elapsed', time.time()-t0)
print(r.stdout[-4000:])
print('STDERR', r.stderr[-1500:])
