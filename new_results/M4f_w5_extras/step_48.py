WORK='<project>/research/08_theory/M4f_w5_extras'
src = r'''"""M4f step 3g: final mode-A polish. Shakes (3-4 coords), finer fractional steps, dual seeds.
"""
import pickle, sys, random, time
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
sol = NS.solve(list(NS.S0), 'opt6')
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
def stats(v):
    v = reduced(v)
    nz = [a for a in v if a]
    terms = sum(1 for p,q in PAIRS if v[p] or v[q])
    return (max(abs(a.numerator) for a in nz), lcm(*[a.denominator for a in nz]), terms)
score = stats

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
        red = flint.fmpz_mat(W).lll(gram='exact')
        self.B = [[int(red[i,j]) for j in range(NATOM)] for i in range(self.k)]
    def t_to_v(self, tt):
        v = [F(w0s, self.D) for w0s in self.w0]
        for i,ti in enumerate(tt):
            if ti:
                for s in range(NATOM):
                    if self.B[i][s]: v[s] += ti*F(self.B[i][s], self.D)
        return v
    def coords_of(self, v):
        diff = [a - F(w0s, self.D) for a, w0s in zip(v, self.w0)]
        M = [[str(F(self.B[j][s], self.D)) for j in range(self.k)] + [str(diff[s])] for s in range(NATOM)]
        mat, rank = flint.fmpq_mat(M).rref()
        pivots = [next((j for j in range(mat.ncols()) if mat[i,j]), None) for i in range(rank)]
        tt = [F(0)]*self.k
        for i,p in enumerate(pivots):
            if p is not None and p < self.k: tt[p] = F(str(mat[i,self.k]))
        return tt

def descend(sp, tt, passes=25,
            steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2),F(1,3),F(-1,3),F(2,3),F(-2,3),
                   F(1,6),F(-1,6),F(1,7),F(-1,7),F(1,12),F(-1,12),F(5,6),F(-5,6))):
    best_v = sp.t_to_v(tt); best_s = score(best_v); tt = list(tt)
    for _ in range(passes):
        improved = False
        for i in range(sp.k):
            for st in steps:
                t2 = list(tt); t2[i] += st
                v2 = sp.t_to_v(t2); s2 = score(v2)
                if s2 < best_s: best_s, tt, best_v = s2, t2, v2; improved = True
        for i in range(sp.k):
            for j in range(i+1, sp.k):
                for si in (F(1),F(-1),F(1,2),F(-1,2)):
                    for sj in (F(1),F(-1),F(1,2),F(-1,2)):
                        t2 = list(tt); t2[i] += si; t2[j] += sj
                        v2 = sp.t_to_v(t2); s2 = score(v2)
                        if s2 < best_s: best_s, tt, best_v = s2, t2, v2; improved = True
        if not improved: break
    return best_s, best_v, tt

def main(budget=1150):
    results = pickle.loads((HERE/'subset_exact.pkl').read_bytes())
    Sfull = (0,1,2,3,4,5,6)
    r = [r for r in results if r['S']==Sfull][0]
    sp = Space(r['v0'], r['NSimg'])
    out5 = pickle.loads((HERE/'pareto.pkl').read_bytes())
    seeds = [sp.coords_of(pickle.loads((HERE/'state_opt.pkl').read_bytes())['v']),
             sp.coords_of(out5['B'][1])]
    best = None
    for sd in seeds:
        sc, bv, tt = descend(sp, sd)
        if best is None or sc < best[0]: best = (sc, bv, tt)
    log(f'seed-descended: {best[0]}')
    rng = random.Random(20250601)
    t0 = time.monotonic(); since = 0; rs = 0
    while time.monotonic()-t0 < budget and since < 120:
        rs += 1
        t2 = list(best[2])
        for _ in range(rng.choice([1,2,2,3,4])):
            i = rng.randrange(sp.k)
            t2[i] += F(rng.choice([1,-1,2,-2,3,-3]))*rng.choice([F(1),F(1),F(1,2),F(1,3),F(1,6)])*rng.choice([1,1,1,2])
        sc2, bv2, tt2 = descend(sp, t2, passes=10,
                                steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2),F(1,3),F(-1,3),F(2,3),F(-2,3),F(1,6),F(-1,6),F(1,7),F(-1,7)))
        if sc2 < best[0]:
            best = (sc2, bv2, tt2); since = 0
            log(f'restart {rs}: {sc2}  ({time.monotonic()-t0:.0f}s)')
            pickle.dump(dict(score=sc2, v=bv2, tt=[str(x) for x in tt2]), open(HERE/'state_opt2.pkl','wb'))
        else: since += 1
    log(f'FINAL: {best[0]} after {rs} restarts')
    pickle.dump(dict(score=best[0], v=best[1], tt=[str(x) for x in best[2]]), open(HERE/'state_opt2.pkl','wb'))

if __name__=='__main__':
    main(int(sys.argv[1]) if len(sys.argv)>1 else 1150)
'''
open(f'{WORK}/optimize6.py','w').write(src)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','optimize6.py','1000'],
                   cwd=WORK, capture_output=True, text=True, timeout=1300)
print('elapsed', time.time()-t0)
print(r.stdout[-2500:])
print('STDERR', r.stderr[-800:])
