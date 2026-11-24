import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
verifier = '''"""verify_oct5d3.py — OCT5D3 verifier (< 5 min at >= 300 digits).

Verifies the closed form of the regular octagon's P^(5)(1):
  A. (PROVED, exact symbol algebra) the weight-5 symbol S(W5oct): obstruction
     (-1485+1035sqrt2, 900-540sqrt2); C's pair equals it; T = S(W5)-S(C) lies exactly
     in the span of the 52 symbol-visible M4f atoms (rank 30/32, residual tensor 0).
  B. (NUMERICAL) the true weight parts Wt0..Wt5 of P^(5)(1)/u at >= 300 digits:
     Wt0 = 0 (exact), Wt1 (also exact), Wt2, Wt3, Wt4 closed forms, and
     Wt5 = sum beta_j (52 classical atoms) + sum gamma (9 Li32) + sum (21 planner
     invisible products)  [the Task-3 closure, two-precision LLL, accepted].
  C. (NUMERICAL) assembly: closed form reproduces P^(5)(1)/u and P^(5)(1) (vs p5_800).
  D. atom checks: Li32 convention (OCT5D2's stored values were off by 1/x), duplications,
     K32 value, L(chi8,5).
  E. prints the m=4 local law and the zeta(5), zeta(2)zeta(3) coefficients.
Labels: PROVED = exact arithmetic; NUMERICAL = high-precision. Exit nonzero on failure.
"""
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

DPS = 360

# ================= A. exact symbol algebra =================
print("== A. exact symbol algebra (PROVED) ==")
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
Am = [row[:] + [bvec[i]] for i, row in enumerate(M)]
pivots = []; r = 0
for c in range(nc):
    piv = next((i for i in range(r, nr) if Am[i][c] != ZERO), None)
    if piv is None: continue
    Am[r], Am[piv] = Am[piv], Am[r]
    pv = Am[r][c]
    Am[r] = [x / pv for x in Am[r]]
    for i in range(nr):
        if i != r and Am[i][c] != ZERO:
            f = Am[i][c]
            Am[i] = [Am[i][j] - f * Am[r][j] for j in range(nc + 1)]
    pivots.append(c); r += 1
    if r == nr: break
incons = [i for i in range(nr) if all(Am[i][c] == ZERO for c in range(nc)) and Am[i][nc] != ZERO]
check("T in span of the 52 visible M4f atoms (rank 30 of 32)", not incons and r == 30,
      "rank %d, kernel %d" % (r, nc - r))
x0 = [ZERO] * nc
for i, c in enumerate(pivots): x0[c] = Am[i][nc]
resid = [sum((M[i][c] * x0[c] for c in range(nc)), ZERO) - bvec[i] for i in range(nr)]
check("residual tensor of the symbol fit is exactly zero", all(v == ZERO for v in resid))
nz = [(vis[i][0], x0[i]) for i in range(nc) if x0[i] != ZERO]
print("   classical representative (M4f basis, %d terms), first 3:" % len(nz))
for lab, v in nz[:3]: print("      %-24s %s" % (lab, v))
import sympy as sp
LI32A = [('(b,1)', (B, ONE)), ('(1/sqrt2,1)', (H, ONE)), ('(-b,1)', (-B, ONE)), ('(b,-1)', (B, -ONE)),
         ('(-b,-1)', (-B, -ONE)), ('(1/sqrt2,-1)', (H, -ONE)), ('(1-1/sqrt2,1)', (ONE - H, ONE)),
         ('(2-sqrt2,1)', (2 - R2, ONE)), ('(-b^2,1)', (-B * B, ONE))]
syms = [li32_symbol(x, y) for _, (x, y) in LI32A]
ws = sorted(set().union(*[set(s) for s in syms]))
Ms = sp.Matrix([[sp.Rational(s.get(w, ZERO).p.numerator, s.get(w, ZERO).p.denominator)
                 for s in syms] for w in ws])
check("Li_{3,2} symbol kernel has dim 1 (the K32 direction)", len(Ms.nullspace()) == 1)

# ================= B. numeric weight parts =================
print("== B. weight parts Wt0..Wt5 at %d digits (NUMERICAL) ==" % DPS)
import gengine2p
G1, letters = gengine2p.gtab(DPS, 1250, cachedir=HERE, procs=min(8, os.cpu_count()))
mp.mp.dps = DPS
L0, LP = letters['L0'], letters['LP']
check("G-table self-check G(0^4,p;1) = -Li5(1/p)",
      abs(G1[(L0,) * 4 + (LP,)] + mp.polylog(5, 1 / LP)) < mp.mpf('1e-330'))
r2 = mp.sqrt(2); L2C = mp.log(2); LU = mp.log(1 + r2)
L1, LQ = letters['L1'], letters['LQ']
def q2f(z): return mp.mpf(str(z.p)) + mp.mpf(str(z.q)) * r2
def w2mpf(w): return tuple(L0 if x == ZERO else L1 if x == ONE else LP if x == P_LET else LQ for x in w)
Gv = {}
def gval(w):
    if w not in Gv: Gv[w] = G1[w2mpf(w)]
    return Gv[w]
allrw = dict(rewrite)
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
check("Wt0 is exactly 0 (QS2 arithmetic)", Wt0_exact * 960 == Q(0))
p5 = pickle.load(open(OCT5D + "p5_800.pkl", "rb"))
err = abs(sum(Wt.values()) - mp.mpf(p5['P5over_u']))
check("sum(Wt_k) = P^(5)(1)/u (vs 800-digit input)", err < mp.mpf('1e-330'),
      "residual " + mp.nstr(err, 3))

from oct5d3_atoms import base_atoms, prod, li32_xy, INVISIBLE5
Aatoms = base_atoms()
import coeff_tables as CT
def spec_of(label):
    parts = label.split('*'); out = []
    for p in parts:
        if p.startswith('lam^'): out.append(('lam', int(p[4:])))
        elif p == 'lam': out.append(('lam', 1))
        elif p.startswith('L^'): out.append(('L', int(p[2:])))
        elif p == 'L': out.append(('L', 1))
        elif p == 'pi2': out.append(('pi2', 1))
        elif p == 'pi2^2': out.append(('pi2', 2))
        elif p == 'z3': out.append(('z3', 1))
        elif p.startswith('Li2(b)^'): out.append(('Li2(b)', int(p[7:])))
        else: out.append((p, 1))
    return out
def cf2(a, b):
    fa, fb = F(a), F(b)
    return mp.mpf(fa.numerator) / fa.denominator + (mp.mpf(fb.numerator) / fb.denominator) * r2
for lvl, tbl in [(1, CT.Wt1), (2, CT.Wt2), (3, CT.Wt3), (4, CT.Wt4)]:
    tot = sum(cf2(a, b) * prod(Aatoms, spec_of(lab)) for lab, (a, b) in tbl.items())
    err = abs(tot - Wt[lvl])
    check("Wt%d closed form" % lvl, err < mp.mpf('1e-330'), "residual " + mp.nstr(err, 3))

# Wt5 primitive expansion
from oct5d3_basis import atom_specs
spec_map = {lab: spec for lab, kind, spec in atom_specs()}
LI32V = {}
bb = r2 - 1
for lab, (x, y) in {'Li32(b,1)': (bb, 1), 'Li32(1/sqrt2,1)': (1 / r2, 1), 'Li32(-b,1)': (-bb, 1),
                    'Li32(b,-1)': (bb, -1), 'Li32(-b,-1)': (-bb, -1), 'Li32(1/sqrt2,-1)': (1 / r2, -1),
                    'Li32(1-1/sqrt2,1)': (1 - 1 / r2, 1), 'Li32(2-sqrt2,1)': (2 - r2, 1),
                    'Li32(-b^2,1)': (-bb * bb, 1)}.items():
    LI32V[lab] = li32_xy(x, y)
tot5 = sum(cf2(a, b) * prod(Aatoms, spec_map[lab]) for lab, (a, b) in CT.WT5_BETA.items())
tot5 += sum(cf2(a, b) * LI32V[lab] for lab, (a, b) in CT.WT5_LI32.items())
planner_specs = dict(INVISIBLE5)
tot5 += sum(cf2(a, b) * prod(Aatoms, planner_specs[nm]) for nm, (a, b) in CT.WT5_PLANNER.items())
err = abs(tot5 - Wt[5])
check("Wt5 = classical + Li32 + planner-invisible closure", err < mp.mpf('1e-330'),
      "residual " + mp.nstr(err, 3))

# ================= C. assembly =================
print("== C. assembly of P^(5)(1) ==")
total_over_u = sum((sum(cf2(a, b) * prod(Aatoms, spec_of(lab)) for lab, (a, b) in CT.__dict__[f'Wt{l}'].items())
                    for l in (1, 2, 3, 4)), mp.mpf(0)) + tot5
u = mp.sqrt(2 + r2)
P5 = u * total_over_u
err1 = abs(total_over_u - mp.mpf(p5['P5over_u']))
err2 = abs(P5 - mp.mpf(p5['P5']))
check("closed form reproduces P^(5)(1)/u", err1 < mp.mpf('1e-330'), "residual " + mp.nstr(err1, 3))
check("closed form reproduces P^(5)(1)", err2 < mp.mpf('1e-330'), "residual " + mp.nstr(err2, 3))

# ================= D. atom checks =================
print("== D. atom checks ==")
av = pickle.load(open(OCT5D2 + "atom_vals430.pkl", "rb"))
check("Li32 convention: Li32(b,1) = b * (OCT5D2 stored)",
      abs(Aatoms['Li32(b,1)'] - bb * mp.mpf(av['L32_b1'])) < mp.mpf('1e-300'),
      "residual " + mp.nstr(abs(Aatoms['Li32(b,1)'] - bb * mp.mpf(av['L32_b1'])), 3))
check("duplication Li5(b^4) = 256(Li5(b)+Li5(-b)) + 16 Li5(-b^2)",
      abs(Aatoms['Li5(b^4)'] - 256 * (Aatoms['Li5(b)'] + Aatoms['Li5(-b)']) - 16 * Aatoms['Li5(-b^2)']) < mp.mpf('1e-330'),
      "residual " + mp.nstr(abs(Aatoms['Li5(b^4)'] - 256 * (Aatoms['Li5(b)'] + Aatoms['Li5(-b)']) - 16 * Aatoms['Li5(-b^2)']), 3))
from oct5d3_atoms import k32_value
K32 = k32_value()
check("K32 = Li32(-b^2,1)-8Li32(b,1)-32Li32(-b,1)-16Li32(b,-1)-16Li32(-b,-1) value",
      abs(K32 - mp.mpf('-0.09874710398912142143877262163927535713803')) < mp.mpf('1e-40'),
      mp.nstr(K32, 25))
import json
l5 = json.load(open(ROOT + "research/08_theory/M4f_w5_extras/Lchi8_5.json"))
check("L(chi8,5) vs M4e certified value", abs(Aatoms['Lchi5'] - mp.mpf(l5['value'])) < mp.mpf('1e-300'),
      "residual " + mp.nstr(abs(Aatoms['Lchi5'] - mp.mpf(l5['value'])), 3))

# ================= E. local law + zeta coefficients =================
print("== E. local law m=4 and zeta coefficients (basis-dependent) ==")
s2q = (2 - Q(0, 1)) / 4; s4 = s2q * s2q; s6 = s4 * s2q
A5 = 65 * s6 - 63 * s4 + 17 * s2q - 1
A23 = 45 * s6 - 27 * s4 + 6 * s2q - 1
C5 = 2880 * 4 * s4 * A5; C23 = -960 * 4 * s4 * A23
c = u / 2
print("   C_4^(5)   = (%s)*cos(pi/8) = %s" % (C5, mp.nstr(q2f(C5) * c, 25)))
print("   C_4^(2,3) = (%s)*cos(pi/8) = %s" % (C23, mp.nstr(q2f(C23) * c, 25)))
z5c = CT.WT5_PLANNER['z5']
z23 = tuple(F(x) for x in CT.WT5_PLANNER['z2*z3'])
z23b = tuple(F(x) for x in CT.WT5_PLANNER['z3*pi2'])
z5v = cf2(*z5c); z23v = cf2(z23[0], z23[1]) + 6 * cf2(z23b[0], z23b[1])
print("   zeta(5)        coefficient in P^(5)(1) = u*(%s+%s sqrt2) = %s" % (z5c[0], z5c[1], mp.nstr(u * z5v, 20)))
print("   zeta(2)zeta(3) coefficient (with pi2*z3 folded) = %s" % mp.nstr(u * z23v, 20))
print("   (non-local parts = these minus the local law; basis-dependent)")

print("\ntotal time %.1fs" % (time.time() - T0))
if failures:
    print("FAILURES:", failures); sys.exit(1)
print("ALL CHECKS PASS")
'''
open("verify_oct5d3.py","w").write(verifier)
print("verify_oct5d3.py written, %d bytes" % len(verifier))
