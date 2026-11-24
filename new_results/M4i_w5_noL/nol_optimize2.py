"""M4i step C2: heavy optimisation of the L-free affine space (optimize3-style).

Iteratively reweighted LLL-CVP (Babai nearest plane, weights (1+(|num_s|/maxnum)^4))
+ rational descent + random restarts. Score = (max|num|, lcm, terms) after the
three reductions. Resumable: seeds from state_noL.pkl; saves on improvement.
"""
import json, pickle, sys, random, time
from fractions import Fraction as F
from math import lcm
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
import nice_search as NS

mp.mp.dps = 200
def log(s): print(s, flush=True)

# --- rebuild the L-free space exactly as nol_build.py ---
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
    else: kidx.append(72+int(nm[1:]))
Bcols = [[F(int(i == j)) for i in range(86)] for j in kidx]
M = [[Bcols[j][i] for j in range(72)] + [deps86[d][i] for d in range(14)] for i in range(86)]
rhs = [[F(int(i == j)) for j in range(86)] for i in range(86)]
mat, rank = flint.fmpq_mat([[str(c) for c in M[i]+rhs[i]] for i in range(86)]).rref()
inv = [[F(str(mat[i, 86+j])) for j in range(86)] for i in range(86)]
coords = [[inv[i][c] for i in range(72)] for c in range(86)]
tbar = [F(-relv[i], c0) for i in range(72)]
sol = NS.solve(list(NS.S0), 'nol2')
names = sol['names']; xp0, xq0 = sol['xp'], sol['xq']
K = json.loads((HERE/'nice_ksat_s0.json').read_text())['K']
pool = json.loads((HERE/'pool.json').read_text())['pool']
s0 = set(NS.S0)
extras = [e for e in pool if NS.W.Q(F(e['p']), F(e['q'])) not in s0]
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
allow = list(range(70)) + list(range(72, 86))
ncol = len(allow)
Mx = [[str(coords[c][i]) for c in allow] + [str(tbar[i])] for i in range(72)]
mat, rank = flint.fmpq_mat(Mx).rref()
pivots = [next((j for j in range(mat.ncols()) if mat[i, j]), None) for i in range(rank)]
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

class Space:
    def __init__(self, v0, NSimg):
        self.k = len(NSimg)
        D = 1
        for a in v0: D = lcm(D, a.denominator)
        for nv in NSimg:
            for a in nv: D = lcm(D, a.denominator)
        self.D = D
        self.w0 = [int(a*D) for a in v0]
        self.W = [[int(a*D) for a in nv] for nv in NSimg]
    def lattice(self):
        red = flint.fmpz_mat(self.W).lll(gram='exact')
        return [[int(red[i, j]) for j in range(NATOM)] for i in range(self.k)]
    def t_to_v(self, tt, B):
        v = [F(w0s, self.D) for w0s in self.w0]
        for i, ti in enumerate(tt):
            if ti:
                for s in range(NATOM):
                    if B[i][s]: v[s] += ti*F(B[i][s], self.D)
        return v

def coords_of(sp, B, v):
    diff = [a - F(w0s, sp.D) for a, w0s in zip(v, sp.w0)]
    Mx2 = [[str(F(B[j][s], sp.D)) for j in range(sp.k)] + [str(diff[s])] for s in range(NATOM)]
    m2, r2 = flint.fmpq_mat(Mx2).rref()
    pv = [next((j for j in range(m2.ncols()) if m2[i, j]), None) for i in range(r2)]
    tt = [F(0)]*sp.k
    for i, p in enumerate(pv):
        if p is not None and p < sp.k: tt[p] = F(str(m2[i, sp.k]))
    return tt

def babai(w0, B):
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
            steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2),F(1,3),F(-1,3),F(2,3),F(-2,3),
                   F(1,6),F(-1,6),F(1,7),F(-1,7),F(1,12),F(-1,12),F(5,6),F(-5,6))):
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
                for si in (F(1),F(-1),F(1,2),F(-1,2)):
                    for sj in (F(1),F(-1),F(1,2),F(-1,2)):
                        t2 = list(tt); t2[i] += si; t2[j] += sj
                        v2 = sp.t_to_v(t2, B); s2 = score(v2)
                        if s2 < best_s: best_s, tt, best_v = s2, t2, v2; improved = True
        if not improved: break
    return best_s, best_v, tt

sp = Space(v0, NSimg)
B = sp.lattice()
log(f'space k={sp.k} D={sp.D}')
statef = HERE/'state_noL.pkl'
cur = pickle.loads(statef.read_bytes())
best_sc, best_v = tuple(cur['score']), cur['v']
best_tt = coords_of(sp, B, best_v)
log(f'start: {best_sc}')

# iteratively reweighted CVP rounds
for rnd in range(8):
    v = sp.t_to_v(best_tt, B)
    sc = score(v)
    weights = [1.0 + (abs(float(a))*sp.D/sc[0])**4 for a in v]
    Ww = [[int(round(x*w)) for x, w in zip(row, weights)] for row in sp.W]
    red = flint.fmpz_mat(Ww).lll(gram='exact')
    Bw = [[int(red[i, j]) for j in range(red.ncols())] for i in range(red.nrows())]
    w0w = [int(round(x*w)) for x, w in zip(sp.w0, weights)]
    ttw = babai(w0w, Bw)
    v2 = [F(w0w[s]) for s in range(NATOM)]
    for i, ti in enumerate(ttw):
        for s in range(NATOM):
            if Bw[i][s]: v2[s] -= ti*Bw[i][s]
    v2 = [F(v2[s])/F(int(round(weights[s])))/sp.D for s in range(NATOM)]
    tt2 = coords_of(sp, B, v2)
    sc2, bv2, tt2 = descend(sp, B, tt2, passes=15)
    log(f'reweighted round {rnd}: {sc2}')
    if sc2 < score(sp.t_to_v(best_tt, B)):
        best_tt = tt2; best_v = bv2
sc, bv, tt = descend(sp, B, best_tt)
log(f'after reweight+descend: {sc}')
if sc < tuple(best_sc):
    pickle.dump(dict(score=sc, v=bv, tt=[str(x) for x in tt]), open(statef, 'wb'))
    log('saved improved state_noL.pkl')
best = (sc, bv, tt)
rng = random.Random(31337)
t0 = time.monotonic(); since = 0; rs = 0
budget = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
while time.monotonic()-t0 < budget and since < 1000000:
    rs += 1
    mode = rng.random()
    t2 = list(best[2])
    if mode < 0.6:
        for _ in range(rng.choice([1, 2, 2, 3, 4])):
            i = rng.randrange(sp.k)
            t2[i] += F(rng.choice([1, -1, 2, -2, 3, -3]))*rng.choice([F(1), F(1), F(1, 2), F(1, 3), F(1, 6)])*rng.choice([1, 1, 1, 2])
    else:  # full-vector shake (optimize3 style)
        t2 = [x + F(rng.choice([1, -1, 2, -2]))*rng.choice([F(1), F(1, 2), F(1, 3)]) for x in t2]
    sc2, bv2, tt2 = descend(sp, B, t2, passes=10,
                            steps=(F(1),F(-1),F(2),F(-2),F(1,2),F(-1,2),F(1,3),F(-1,3),F(2,3),F(-2,3),F(1,6),F(-1,6),F(1,7),F(-1,7)))
    if sc2 < best[0]:
        best = (sc2, bv2, tt2); since = 0
        log(f'restart {rs}: {sc2}  ({time.monotonic()-t0:.0f}s)')
        pickle.dump(dict(score=sc2, v=bv2, tt=[str(x) for x in tt2]), open(statef, 'wb'))
    else: since += 1
log(f'FINAL noL2: {best[0]} after {rs} restarts')
pickle.dump(dict(score=best[0], v=best[1], tt=[str(x) for x in best[2]]), open(statef, 'wb'))
