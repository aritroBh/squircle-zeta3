"""M4i step A: explicit L-relation. L(chi8,5) = sum (p+q*sqrt2)*atom, no L on RHS.

Lambda = rowspan(deps86) (14 exact value relations among 86 columns, M4f overcomplete_1200,
same span reconfirmed at 1500). Columns: 0..34 kv, 35..69 s2kv, 70 L, 71 s2L, 72..78 u0..6, 79..85 s2u0..6.
- particular lam0 in Lambda with lam0[70]=0, lam0[71]=1  (s2L = sqrt2*L via other columns)
- Lambda00 = {lam in Lambda : lam[70]=lam[71]=0}  (12-dim freedom of exact zero relations)
Atom map ACOLS (146 slots) exactly as subset_exact.py. vLsqrt2 = e_131 - ACOLS(lam0) has
value sqrt2*L and zero L-slots. swap((p,q))=(q,p/2) divides values by sqrt2 -> vL has value L.
Minimise (max|num|, lcm, terms) of reduce3(vL) over the 12-dim lattice by
LLL(exact gram)+rational descent+shakes, as optimize6/7. Save L_relation.json (R5_best format).
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

# ---------- exact column machinery (as subset_exact.py) ----------
oc = json.loads((HERE/'overcomplete_1200.json').read_text())
deps86 = [[F(c) for c in dep['vec'][1:]] for dep in oc['deps']]
assert len(deps86) == 14 and all(len(d) == 86 for d in deps86)

# particular lam0: lam0[70]=0, lam0[71]=1 ; basis of Lambda00
rows = [[str(deps86[d][70]) for d in range(14)] + ['0'],
        [str(deps86[d][71]) for d in range(14)] + ['1']]
mat, rank = flint.fmpq_mat(rows).rref()
assert rank == 2, f'L/s2L directions in Lambda have rank {rank} != 2'
pivots = [next((j for j in range(15) if mat[i, j]), None) for i in range(rank)]
assert all(p < 14 for p in pivots)
a0 = [F(0)]*14
for i, p in enumerate(pivots): a0[p] = F(str(mat[i, 14]))
free = [j for j in range(14) if j not in pivots]
nabasis = []
for f in free:
    w = [F(0)]*14; w[f] = F(1)
    for i, p in enumerate(pivots): w[p] = -F(str(mat[i, f]))
    nabasis.append(w)
assert len(nabasis) == 12
lam0 = [sum(a0[d]*deps86[d][c] for d in range(14)) for c in range(86)]
assert lam0[70] == 0 and lam0[71] == 1
L00 = []
for na in nabasis:
    la = [sum(na[d]*deps86[d][c] for d in range(14)) for c in range(86)]
    assert la[70] == 0 and la[71] == 0
    L00.append(la)
log(f'particular lam0 + 12-dim Lambda00 constructed (exact)')

# readable u-form: s2L = -sum_{c != 71} lam0[c] col_c
unames = {72+i: f'u{i}' for i in range(7)} | {79+i: f's2u{i}' for i in range(7)}
usup = [(unames.get(c, f'kv{c}' if c < 35 else f's2kv{c-35}'), -lam0[c])
        for c in range(86) if c != 71 and lam0[c] and (c >= 72 or c in (70,))]
log('s2L = sqrt2*L relation, u/L support: ' + str(usup))
kvsup = [(f'kv{c}' if c < 35 else f's2kv{c-35}', -lam0[c]) for c in range(70) if lam0[c]]
log(f'  kv support ({len(kvsup)} entries): {kvsup}')

# ---------- atom map (identical to subset_exact.py) ----------
sol = NS.solve(list(NS.S0), 'Lrel')
names = sol['names']
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
def aimag(lam):
    v = [F(0)]*NATOM
    for c in range(86):
        if lam[c]:
            A = ACOLS[c]
            for s in range(NATOM):
                if A[s]: v[s] += lam[c]*A[s]
    return v

def swapmap(v):
    w = [F(0)]*NATOM
    for i in range(65):
        w[i] = v[65+i]; w[65+i] = v[i]/2
    w[130] = v[131]; w[131] = v[130]/2
    for k in range(7):
        w[132+k] = v[139+k]; w[139+k] = v[132+k]/2
    return w

# reductions (all three), as export_forms.reduce3
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

# base and freedom in L-form
e131 = [F(0)]*NATOM; e131[131] = F(1)
vLsqrt2 = [e131[s] - aimag(lam0)[s] for s in range(NATOM)]
base = swapmap(vLsqrt2)
assert base[130] == 0 and base[131] == 0
G = [swapmap(aimag(la)) for la in L00]
log(f'raw L-relation score (max|num|, lcm, terms): {score(base)}')

# ---------- optimise over the 12-dim lattice (optimize6/7 machinery) ----------
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
        M = [[str(F(self.B[j][s], self.D)) for j in range(self.k)] + [str(diff[s])] for s in range(NATOM)]
        mat, rank = flint.fmpq_mat(M).rref()
        pivots = [next((j for j in range(mat.ncols()) if mat[i, j]), None) for i in range(rank)]
        tt = [F(0)]*self.k
        for i, p in enumerate(pivots):
            if p is not None and p < self.k: tt[p] = F(str(mat[i, self.k]))
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

sp = Space(base, G)
log(f'space: k={sp.k} D={sp.D}')
sc, bv, tt = descend(sp, [F(0)]*sp.k)
log(f'descended from raw: {sc}')
rng = random.Random(20261022)
t0 = time.monotonic(); since = 0; rs = 0
budget = int(sys.argv[1]) if len(sys.argv) > 1 else 420
best = (sc, bv, tt)
while time.monotonic()-t0 < budget and since < 40:
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
    else: since += 1
log(f'FINAL L-relation score: {best[0]} after {rs} restarts')
vL = best[1]
assert vL[130] == 0 and vL[131] == 0

# ---------- export ----------
extrasQ = [NS.W.Q(F(e['p']), F(e['q'])) for e in extras]
def export(v, path):
    v = reduce3(v)
    terms = {}
    for i, nm in enumerate(names):
        a, b = v[i], v[65+i]
        if a or b: terms[nm] = [str(a), str(b)]
    for k in range(7):
        a, b = v[132+k], v[139+k]
        if a or b: terms[f'Li5({extrasQ[k]})'] = [str(a), str(b)]
    a, b = v[130], v[131]
    if a or b: terms['Lchi8_5'] = [str(a), str(b)]
    assert 'Lchi8_5' not in terms, 'L still present!'
    dens = [F(x).denominator for pq in terms.values() for x in pq]
    nums = [abs(F(x).numerator) for pq in terms.values() for x in pq]
    li5 = [nm for nm in terms if nm.startswith('Li5')]
    log(f'{path}: terms={len(terms)} max|num|={max(nums)} lcm={lcm(*dens)} Li5 args={len(li5)}')
    json.dump({'source': 'M4i L-relation: L(chi8,5) = sum (p+q*sqrt2)*atom, no L on RHS',
               'terms': terms}, open(HERE/path, 'w'), indent=1)
    return terms
terms = export(vL, 'L_relation.json')
pickle.dump(dict(score=best[0], v=vL, tt=[str(x) for x in best[2]]),
            open(HERE/'lrel_state.pkl', 'wb'))
log('saved L_relation.json and lrel_state.pkl')
# 10 largest coefficients
allc = sorted(((max(abs(F(p).numerator), 1), nm, p, q) for nm, (p, q) in terms.items()),
              reverse=True)
log('largest coefficients:')
for mx, nm, p, q in allc[:10]:
    log(f'  |num|={mx}  {nm}: [{p}, {q}]')
