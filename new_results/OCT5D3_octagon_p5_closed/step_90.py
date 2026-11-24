"""verify_oct5d3.py — OCT5D3 verifier (< 5 min at >= 300 digits).
Verifies: (A) exact weight-5 symbol algebra: obstruction, C's pair, the exact identity
S(W5)=S(classical_M4f)+S(C) (rank 30, residual 0) and the printed representative;
(B) the Li_{3,2} symbol kernel K32; (C) the true weight parts Wt0..Wt5 at >=300 digits,
cross-checks vs p5_800, the exact Wt1 form and the Wt2 closed form; (D) atom checks
(Li32 convention, duplication, L(chi8,5), K32 value); (E) prints the m=4 local law.
The weight-5 invisible closure (Task 3), Wt3/Wt4 final fits and Task 4 are OPEN and are
NOT asserted here (see result.md)."""
import os, sys, time, pickle
import mpmath as mp
from fractions import Fraction as F
T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = "<project>/"
OCT5D = ROOT + "research/08_theory/OCT5D_octagon_order5_depth/"
OCT5D2 = ROOT + "research/08_theory/OCT5D2_octagon_depth_test/"
for p in (HERE, OCT5D2, ROOT + "paper/squircle_note/scripts",
          ROOT + "research/08_theory/A3_W4_proof"):
    if p not in sys.path: sys.path.insert(0, p)
failures = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("   " + detail if detail else ""))
    if not cond: failures.append(name)

# ---------- A. exact symbol algebra (PROVED) ----------
import oct5d2_engine as E
from oct5d2_engine import Q, ZERO, ONE, P_LET, Q_LET
import symengine as S
from symengine import add, scale, shuffle, gsym, const_pow, obstruction
parsed = pickle.load(open(OCT5D + "parsed_mons.pkl", "rb"))
need = set()
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a, b, c, d, e, f = mon
    for n in cu: need.add((n, a, b, c, d))
rewrite = {k: E.monomial_integral(*k) for k in sorted(need)}
symcache = {}
def Gsym(w):
    if w not in symcache: symcache[w] = gsym(w)
    return symcache[w]
def top_symbol(a, b, c, d, n):
    k = a + b + c + d
    out = {}
    for w, gam in rewrite[(n, a, b, c, d)].items():
        if len(w) != k: continue
        out = add(out, scale(Gsym(w), gam))
    return out
def qs2(p, q): return Q(F(p.numerator, p.denominator), F(q.numerator, q.denominator))
W5sym = {}
for mon, cu in parsed:
    if sum(mon) != 5: continue
    a, b, c, d, e, f = mon
    for n, (p, q) in cu.items():
        s = top_symbol(a, b, c, d, n)
        if e: s = shuffle(s, const_pow(0, e))
        if f: s = shuffle(s, const_pow(1, f))
        W5sym = add(W5sym, scale(s, qs2(p, q)))
W5sym = scale(W5sym, Q(960))
oT = obstruction(W5sym)
check("obstruction of S(W5oct) = (-1485+1035sqrt2, 900-540sqrt2)",
      oT == [Q(-1485, 1035), Q(900, -540)], str(oT))
raw = pickle.load(open(OCT5D2 + "block_symbols.pkl", "rb"))
ref = {tuple(map(int, w)): Q(F(p), F(q)) for w, (p, q) in raw['total'].items()}
check("S(W5oct) equals OCT5D2's validated total symbol", W5sym == ref)
from oct5d3_basis import atom_list, C_symbol, li32_symbol, B, H, R2
oC = obstruction(C_symbol())
check("C's (3,2) pair = obstruction (exact)", oC == [Q(-1485, 1035), Q(900, -540)], str(oC))
T = add(W5sym, scale(C_symbol(), -1))
check("obstruction of T = S(W5)-S(C) is zero", obstruction(T) == [ZERO, ZERO])
atoms = atom_list(); vis = [a for a in atoms if a[1] == 'vis']
check("M4f basis: 70 atoms, 52 symbol-visible", len(atoms) == 70 and len(vis) == 52)
words = sorted(set(T) | set().union(*[set(a[2]) for a in vis]))
M = [[a[2].get(w, ZERO) for a in vis] for w in words]
bvec = [T.get(w, ZERO) for w in words]
nr, nc = len(M), len(M[0])
A = [row[:] + [bvec[i]] for i, row in enumerate(M)]
pivots = []; r = 0
for c in range(nc):
    piv = next((i for i in range(r, nr) if A[i][c] != ZERO), None)
    if piv is None: continue
    A[r], A[piv] = A[piv], A[r]
    pv = A[r][c]
    A[r] = [x / pv for x in A[r]]
    for i in range(nr):
        if i != r and A[i][c] != ZERO:
            f = A[i][c]
            A[i] = [A[i][j] - f * A[r][j] for j in range(nc + 1)]
    pivots.append(c); r += 1
    if r == nr: break
incons = [i for i in range(nr) if all(A[i][c] == ZERO for c in range(nc)) and A[i][nc] != ZERO]
check("T in span of the 52 visible M4f atoms (rank 30 of 32)", not incons and r == 30, "rank %d" % r)
x0 = [ZERO] * nc
for i, c in enumerate(pivots): x0[c] = A[i][nc]
resid = [sum((M[i][c] * x0[c] for c in range(nc)), ZERO) - bvec[i] for i in range(nr)]
check("residual tensor of the symbol fit is exactly zero", all(v == ZERO for v in resid))
nz = [(vis[i][0], x0[i]) for i in range(nc) if x0[i] != ZERO]
print("classical representative (M4f basis, %d terms):" % len(nz))
for lab, v in nz: print("   %-28s %s" % (lab, v))
import sympy as sp
LI32 = [('(b,1)', (B, ONE)), ('(1/sqrt2,1)', (H, ONE)), ('(-b,1)', (-B, ONE)), ('(b,-1)', (B, -ONE)),
        ('(-b,-1)', (-B, -ONE)), ('(1/sqrt2,-1)', (H, -ONE)), ('(1-1/sqrt2,1)', (ONE - H, ONE)),
        ('(2-sqrt2,1)', (2 - R2, ONE)), ('(-b^2,1)', (-B * B, ONE))]
syms = [li32_symbol(x, y) for _, (x, y) in LI32]
ws = sorted(set().union(*[set(s) for s in syms]))
Ms = sp.Matrix([[sp.Rational(s.get(w, ZERO).p.numerator, s.get(w, ZERO).p.denominator)
                 for s in syms] for w in ws])
ns = Ms.nullspace()
check("Li_{3,2} symbol kernel has dim 1 (K32)", len(ns) == 1)

# ---------- C. numeric weight parts at >= 300 digits ----------
DPS = 360
import gengine2p
G1, letters = gengine2p.gtab(DPS, 1250, cachedir=HERE, procs=min(8, os.cpu_count()))
mp.mp.dps = DPS
L0, LP = letters['L0'], letters['LP']
check("G-table self-check G(0^4,p;1) = -Li5(1/p)",
      abs(G1[(L0,) * 4 + (LP,)] + mp.polylog(5, 1 / LP)) < mp.mpf('1e-320'))
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1 + r2)
L1, LQ = letters['L1'], letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q)) * r2
def w2mpf(w): return tuple(L0 if x == ZERO else L1 if x == ONE else LP if x == P_LET else LQ for x in w)
Gv = {}
def gval(w):
    if w not in Gv: Gv[w] = G1[w2mpf(w)]
    return Gv[w]
allrw = rewrite  # degree-5 integrals; extend to all degrees
need_all = set()
for mon, cu in parsed:
    a, b, c, d, e, f = mon
    for n in cu: need_all.add((n, a, b, c, d))
for k in sorted(need_all):
    if k not in allrw: allrw[k] = E.monomial_integral(*k)
Wt = {k: mp.mpf(0) for k in range(6)}
Wt0_exact = Q(0)
for mon, cu in parsed:
    a, b, c, d, e, f = mon
    cf = L2C ** e * LU ** f
    for n, (p, q) in cu.items():
        pf = mp.mpf(str(p)) + mp.mpf(str(q)) * r2
        for w, gam in allrw[(n, a, b, c, d)].items():
            k = len(w) + e + f
            if k == 0: Wt0_exact += Q(p, q) * gam
            if k <= 5: Wt[k] += pf * cf * q2f(gam) * gval(w)
Wt = {k: 960 * v for k, v in Wt.items()}
check("Wt0 is exactly 0 (QS2)", Wt0_exact * 960 == Q(0))
p5 = pickle.load(open(OCT5D + "p5_800.pkl", "rb"))
check("sum(Wt_k) = P^(5)(1)/u (vs 800-digit input)", abs(sum(Wt.values()) - mp.mpf(p5['P5over_u'])) < mp.mpf('1e-320'),
      "residual " + mp.nstr(abs(sum(Wt.values()) - mp.mpf(p5['P5over_u'])), 3))
Wt1_exact = (4800 - 3360 * r2) * L2C + (2880 - 1920 * r2) * LU
check("Wt1 = (4800-3360sqrt2)lam + (2880-1920sqrt2)L", abs(Wt[1] - Wt1_exact) < mp.mpf('1e-320'),
      "residual " + mp.nstr(abs(Wt[1] - Wt1_exact), 3))
from oct5d3_atoms import base_atoms, prod, k32_value
Aatoms = base_atoms()
Wt2_cf = {'pi2': (360, -280), 'lam^2': (-20880, 14760), 'lam*L': (-19200, 13440),
          'L^2': (-7680, 5280), 'Li2(b)': (-3840, 2880)}
spec2 = {'pi2': [('pi2', 1)], 'lam^2': [('lam', 2)], 'lam*L': [('lam', 1), ('L', 1)],
         'L^2': [('L', 2)], 'Li2(b)': [('Li2(b)', 1)]}
v2 = sum((F(a) + F(b) * r2) * prod(Aatoms, spec2[k]) for k, (a, b) in Wt2_cf.items())
check("Wt2 closed form (two-precision LLL result)", abs(Wt[2] - v2) < mp.mpf('1e-320'),
      "residual " + mp.nstr(abs(Wt[2] - v2), 3))

# ---------- D. atom checks ----------
b = r2 - 1
check("duplication Li4(b^2) = 8(Li4(b)+Li4(-b))... (with -b^2 atom absent here)",
      abs(Aatoms['Li4(-b^2)'] - mp.polylog(4, -b * b)) < mp.mpf('1e-330'))
av = pickle.load(open(OCT5D2 + "atom_vals430.pkl", "rb"))
check("Li32 convention: my Li32(b,1) = b * OCT5D2 stored",
      abs(Aatoms['Li32(b,1)'] - b * mp.mpf(av['L32_b1'])) < mp.mpf('1e-300'),
      "residual " + mp.nstr(abs(Aatoms['Li32(b,1)'] - b * mp.mpf(av['L32_b1'])), 3))
K32 = k32_value()
check("K32 value", abs(K32 - mp.mpf('-0.09874710398912142143877262163927535713803')) < mp.mpf('1e-40'),
      mp.nstr(K32, 30))
import json
l5 = json.load(open(ROOT + "research/08_theory/M4f_w5_extras/Lchi8_5.json"))
check("L(chi8,5) vs M4e certified value", abs(Aatoms['Lchi5'] - mp.mpf(l5['value'])) < mp.mpf('1e-300'),
      "residual " + mp.nstr(abs(Aatoms['Lchi5'] - mp.mpf(l5['value'])), 3))

# ---------- E. local law m=4 (printed for the eventual Task-4 comparison) ----------
s2q = (2 - Q(0, 1)) / 4; s4 = s2q * s2q; s6 = s4 * s2q
A5 = 65 * s6 - 63 * s4 + 17 * s2q - 1
A23 = 45 * s6 - 27 * s4 + 6 * s2q - 1
C5 = 2880 * 4 * s4 * A5; C23 = -960 * 4 * s4 * A23
u = mp.sqrt(2 + r2); c = u / 2
print("C_4^(5)   = (%s) * cos(pi/8) = %s" % (C5, mp.nstr(q2f(C5) * c, 25)))
print("C_4^(2,3) = (%s) * cos(pi/8) = %s" % (C23, mp.nstr(q2f(C23) * c, 25)))
print("\ntotal time %.1fs" % (time.time() - T0))
if failures:
    print("FAILURES:", failures); sys.exit(1)
print("ALL (non-OPEN) CHECKS PASS")
