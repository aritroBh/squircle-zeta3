import os
dst = '<project>/research/08_theory/M4i_w5_noL'
with open(dst+'/progress.log','a') as fh:
    fh.write('step_19(lrel_verify): symbol=0 PASS; 13 Li5 cobracket-0 PASS; Hurwitz vs file 4.24e-1656; RHS-L(Hurwitz)=2.78e-1741 fresh 1620d.\n')

script = r'''"""M4i step C: W5 WITHOUT L. Build the L-free affine solution space and optimise.

Allowed columns: kv(0..69) + u0..6,s2u0..6 (72..85); NO L/s2L (70,71). 84 columns, rank 72
(drop_L.py confirmed tbar in span), so 12-dim nullspace. v0 = vpart + sum x0 ACOLS.
Optimise (max|num|, lcm, terms) after the three reductions (b^2, -1, -1/sqrt2) with
LLL(exact gram) + rational descent + shake restarts, as optimize6/7. Resumable:
seeds = rref particular, R5_best with L substituted by the M4i L-relation, saved state.
"""
import json, pickle, sys, random, time
from fractions import Fraction as F
from math import lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import nice_search as NS

def log(s): print(s, flush=True)

# ---------- exact coordinate machinery (as subset_exact.py / drop_L.py) ----------
oc = json.loads((HERE/'overcomplete_1200.json').read_text())
deps86 = [[F(c) for c in dep['vec'][1:]] for dep in oc['deps']]
rel = oc['target_rels'][0]; c0 = F(rel['vec'][0]); relv = [F(c) for c in rel['vec'][1:]]
kept = oc['kept']
kidx = []
for nm in kept:
    if nm == 'L': kidx.append(70)
    elif nm == 's2L': kidx.append(71)
    elif nm.startswith('s2kv'): kidx.append(35+int(nm[4:]))
    elif nm.startswith('kv'): kidx.append(int(nm[2:]))
    elif nm.startswith('s2u'): kidx.append(79+int(nm[3:]))
    elif nm.startswith('u'): kidx.append(72+int(nm[1:]))
    else: raise ValueError(nm)
Bcols = [[F(int(i == j)) for i in range(86)] for j in kidx]
M = [[Bcols[j][i] for j in range(72)] + [deps86[d][i] for d in range(14)] for i in range(86)]
rhs = [[F(int(i == j)) for j in range(86)] for i in range(86)]
mat, rank = flint.fmpq_mat([[str(c) for c in M[i]+rhs[i]] for i in range(86)]).rref()
assert rank == 86
inv = [[F(str(mat[i, 86+j])) for j in range(86)] for i in range(86)]
coords = [[inv[i][c] for i in range(72)] for c in range(86)]
tbar = [F(-relv[i], c0) for i in range(72)]

sol = NS.solve(list(NS.S0), 'nol')
names = sol['names']; xp0, xq0 = sol['xp'], sol['xq']
K = json.loads((HERE/'nice_ksat_s0.json').read_text())['K']
pool = json.loads((HERE/'pool.json').read_text())['pool']
s0 = set(NS.S0)
extras = [e for e in pool if NS.W.Q(F(e['p']), F(e['q'])) not in s0]
assert len(extras) == 7
NATOM = 146
def acol(c):
    v = [F(0)]*NATOM
    if c < 35:
        for i, a in enumerate(K[c]): v[i] = F(a)
    elif c < 70:
        for i, a in enumerate(K[c-35]): v[65+i] = F(a)
    elif c == 70: v[130] = F(1)
    elif c == 71: v[131] = F(1)
    elif c < 79:
        i = c-72; e = extras[i]
        v[132+i] = F(1)
        for j, a in enumerate(e['xp']): v[j] -= F(a)
        for j, a in enumerate(e['xq']): v[65+j] -= F(a)
    else:
        i = c-79; e = extras[i]
        v[139+i] = F(1)
        for j, a in enumerate(e['xq']): v[j] -= 2*F(a)
        for j, a in enumerate(e['xp']): v[65+j] -= F(a)
    return v
ACOLS = [acol(c) for c in range(86)]
vpart = [F(0)]*NATOM
for i, a in enumerate(xp0): vpart[i] = a
for i, a in enumerate(xq0): vpart[65+i] = a

# ---------- L-free solve ----------
allow = list(range(70)) + list(range(72, 86))
ncol = len(allow)
Mx = [[str(coords[c][i]) for c in allow] + [str(tbar[i])] for i in range(72)]
mat, rank = flint.fmpq_mat(Mx).rref()
pivots = [next((j for j in range(mat.ncols()) if mat[i, j]), None) for i in range(rank)]
assert all(p is None or p < ncol for p in pivots), 'target NOT in L-free span'
free = [j for j in range(ncol) if j not in pivots]
x0 = [F(0)]*ncol
for i, p in enumerate(pivots):
    if p is not None and p < ncol: x0[p] = F(str(mat[i, ncol]))
ns = []
for f in free:
    w = [F(0)]*ncol; w[f] = F(1)
    for i, p in enumerate(pivots):
        if p is not None and p < ncol: w[p] = -F(str(mat[i, f]))
    ns.append(w)
log(f'L-free space: {ncol} columns, nullity={len(ns)}')
v0 = list(vpart)
for j, c in enumerate(allow):
    if x0[j]:
        A = ACOLS[c]
        for s in range(NATOM):
            if A[s]: v0[s] += x0[j]*A[s]
NSimg = []
for w in ns:
    vi = [F(0)]*NATOM
    for j, c in enumerate(allow):
        if w[j]:
            A = ACOLS[c]
            for s in range(NATOM):
                if A[s]: vi[s] += w[j]*A[s]
    NSimg.append(vi)
assert v0[130] == 0 and v0[131] == 0
assert all(nv[130] == 0 and nv[131] == 0 for nv in NSimg)

# ---------- reductions + score (three reductions, as briefed) ----------
i_b2 = names.index('Li5((3-2*sqrt2))'); i_bp = names.index('Li5((-1+1*sqrt2))')
i_bm = names.index('Li5((1-1*sqrt2))'); i_m1 = names.index('Li5(-1)')
i_z5 = names.index('z5'); i_r2m = names.index('Li5(-1/2*sqrt2)')
i_r2p = names.index('Li5(1/2*sqrt2)'); i_half = names.index('Li5(1/2)')
def reduce3(v):
    v = list(v)
    for off in (0, 65):
        a = v[i_r2m+off]
        if a: v[i_half+off] += F(a, 16); v[i_r2p+off] -= a; v[i_r2m+off] = F(0)
        a = v[i_b2+off]
        if a: v[i_bp+off] += 16*a; v[i_bm+off] += 16*a; v[i_b2+off] = F(0)
        a = v[i_m1+off]
        if a: v[i_z5+off] += F(-15, 16)*a; v[i_m1+off] = F(0)
    return v
PAIRS = [(i, i+65) for i in range(65)] + [(130, 131)] + [(132+i, 139+i) for i in range(7)]
def stats(v):
    v = reduce3(v)
    nz = [a for a in v if a]
    terms = sum(1 for p, q in PAIRS if v[p] or v[q])
    return (max(abs(a.numerator) for a in nz), lcm(*[a.denominator for a in nz]), terms)
score = stats
log(f'rref particular score: {score(v0)}')

# ---------- seed: R5_best with L substituted by the M4i L-relation ----------
vb = pickle.loads((HERE/'state_opt2.pkl').read_bytes())['v']
vL = pickle.loads((HERE/'lrel_state.pkl').read_bytes())['v']
def swapmap(v):
    w = [F(0)]*NATOM
    for i in range(65):
        w[i] = v[65+i]; w[65+i] = v[i]/2
    w[130] = v[131]; w[131] = v[130]/2
    for k in range(7):
        w[132+k] = v[139+k]; w[139+k] = v[132+k]/2
    return w
pL, qL = vb[130], vb[131]
vswap = swapmap(vL)   # value L/sqrt2 -> 2*vswap has value sqrt2 L
v_sub = [vb[s] + pL*vL[s] + 2*qL*vswap[s] for s in range(NATOM)]
v_sub[130] = F(0); v_sub[131] = F(0)
log(f'substituted R5_best score: {score(v_sub)}')

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
        self.B = [[int(red[i, j]) for j in range(NATOM)] for i in range(self.k)]
    def t_to_v(self, tt):
        v = [F(w0s, self.D) for w0s in self.w0]
        for i, ti in enumerate(tt):
            if ti:
                for s in range(NATOM):
                    if self.B[i][s]: v[s] += ti*F(self.B[i][s], self.D)
        return v
    def coords_of(self, v):
        diff = [a - F(w0s, self.D) for a, w0s in zip(v, self.w0)]
        Mx2 = [[str(F(self.B[j][s], self.D)) for j in range(self.k)] + [str(diff[s])] for s in range(NATOM)]
        m2, r2 = flint.fmpq_mat(Mx2).rref()
        pv = [next((j for j in range(m2.ncols()) if m2[i, j]), None) for i in range(r2)]
        tt = [F(0)]*self.k
        for i, p in enumerate(pv):
            if p is not None and p < self.k: tt[p] = F(str(m2[i, self.k]))
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

sp = Space(v0, NSimg)
log(f'space: k={sp.k} D={sp.D}')
seeds = [[F(0)]*sp.k, sp.coords_of(v_sub)]
statef = HERE/'state_noL.pkl'
if statef.exists():
    seeds.append(sp.coords_of(pickle.loads(statef.read_bytes())['v']))
best = None
for sd in seeds:
    sc, bvv, tt = descend(sp, sd)
    log(f'seed descended: {sc}')
    if best is None or sc < best[0]: best = (sc, bvv, tt)
rng = random.Random(20251023)
t0 = time.monotonic(); since = 0; rs = 0
budget = int(sys.argv[1]) if len(sys.argv) > 1 else 1200
while time.monotonic()-t0 < budget and since < 120:
    rs += 1
    t2 = list(best[2])
    for _ in range(rng.choice([1, 2, 2, 3, 4])):
        i = rng.randrange(sp.k)
        t2[i] += F(rng.choice([1, -1, 2, -2, 3, -3]))*rng.choice([F(1), F(1), F(1, 2), F(1, 3), F(1, 6)])*rng.choice([1, 1, 1, 2])
    sc2, bv2, tt2 = descend(sp, t2, passes=10,
                            steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2),F(1,3),F(-1,3),F(2,3),F(-2,3),F(1,6),F(-1,6),F(1,7),F(-1,7)))
    if sc2 < best[0]:
        best = (sc2, bv2, tt2); since = 0
        log(f'restart {rs}: {sc2}  ({time.monotonic()-t0:.0f}s)')
        pickle.dump(dict(score=sc2, v=bv2, tt=[str(x) for x in tt2]), open(statef, 'wb'))
    else: since += 1
log(f'FINAL noL: {best[0]} after {rs} restarts')
pickle.dump(dict(score=best[0], v=best[1], tt=[str(x) for x in best[2]]), open(statef, 'wb'))
'''
open(dst+'/nol_build.py','w').write(script)
print('written nol_build.py', len(script))
