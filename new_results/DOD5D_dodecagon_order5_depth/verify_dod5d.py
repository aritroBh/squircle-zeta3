"""verify_dod5d.py — DOD5D depth-test verifier (< 5 min).

Reproduces/validates, from the read-only inputs and the work-dir deliverables:
  (A) value: engine P^(5)(1) at dps 60 vs the stored 170-digit value; P(1) = 12 sqrt(26-15 r3);
      P‴(1) vs DOD1_dodecagon/P3_630.pkl.
  (B) exact structure: re-parse of E5 (read-only step_02_full.pkl) + exact Q(sqrt3) expansion,
      compared EXACTLY to E5_dod_expanded.pkl; spot check vs symeng.E5_sym; itab140 spot check
      (16 monomials fresh quadrature) + full itab assembly vs engine (>= 60 digits).
  (C) rewrite + G-table: 60 random monomial rewrites recomputed exactly; G-table self-checks
      (quadrature, shuffle, log, polylog); full gtab assembly vs engine (>= 60 digits).
  (D) symbols: 40 random g_symbols recomputed; full total symbol re-assembled (grouped fast path)
      and compared EXACTLY to dod5d_symbols.pkl; obstruction: total table (24 pairings, nonzero),
      local = 0, block sum = total, mirror equalities; controls: Li_5 at 13 S-units -> 0,
      300 shuffle products -> 0, M4c W5 -> (84, -38 sqrt2), octagon -> (-1485+1035 sqrt2,
      900-540 sqrt2); G-symbol identities; 24-term Li_3,2/Li_4,1 correction residual exactly 0.
Exits nonzero on any failure.
"""
import os, sys, time, pickle, json, random
from fractions import Fraction as F
T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
THEORY = "<project>/research/08_theory/"
for p in (HERE, THEORY+"E5c_nonlocal_z5", "<project>/paper/squircle_note/scripts"):
    if p not in sys.path: sys.path.insert(0, p)
failures = []
def check(name, cond, detail=""):
    print(("PASS  " if cond else "FAIL  ") + name + ("   " + str(detail) if detail else ""), flush=True)
    if not cond: failures.append(name)

# ---------------- (A) value ----------------
import mpmath as mp
import dod_engine
res70, diff70 = dod_engine.run(70, None, logf=lambda s: None)
mp.mp.dps = 190
P5_170 = mp.mpf("-62.53724775498639535309129282205803572684769277410373262688145818074695731422874528214220721559218128482460652943087962941126065293896549327018536685224048730049659698854")
check("A1 value: dps-70 engine P^(5)(1) matches 170-digit value",
      abs(mp.mpf(str(res70[5])) - P5_170) < mp.mpf('1e-64'), mp.nstr(abs(mp.mpf(str(res70[5]))-P5_170), 3))
P3_630 = pickle.load(open(THEORY+"DOD1_dodecagon/P3_630.pkl","rb"))
check("A2 P‴(1) vs P3_630.pkl", abs(mp.mpf(str(res70[3])) - P3_630['P3']) < mp.mpf('1e-64'),
      mp.nstr(abs(mp.mpf(str(res70[3])) - P3_630['P3']), 3))
check("A3 P(1) = 12 sqrt(26-15 sqrt3)", abs(mp.mpf(str(res70[0])) - 12*mp.sqrt(26-15*mp.sqrt(3))) < mp.mpf('1e-64'),
      mp.nstr(abs(mp.mpf(str(res70[0])) - 12*mp.sqrt(26-15*mp.sqrt(3))), 3))

# ---------------- (B) exact structure ----------------
import sympy as sp
_t_, _al, _be, _Ls, _hs, _s0s = sp.symbols('t alpha beta L h s0', real=True)
_l0s, _l1s = sp.symbols('l0 l1')
_taus = sp.symbols('tau1:6'); _gams = sp.symbols('gam1:6')
_loc = {'t':_t_,'alpha':_al,'beta':_be,'L':_Ls,'h':_hs,'s0':_s0s,'l0':_l0s,'l1':_l1s}
for _i,_x in enumerate(_taus): _loc[f'tau{_i+1}']=_x
for _i,_x in enumerate(_gams): _loc[f'gam{_i+1}']=_x
_D = pickle.load(open(THEORY+"E5_edge_law/steps/step_02_full.pkl",'rb'))
E5full = sp.sympify(_D['E5'], locals=_loc)
allvars = [_t_,_al,_be,_Ls,_hs,_s0s,_l0s,_l1s]+list(_taus)+list(_gams)
E5exp = sp.expand(E5full)
terms = []
for term in sp.Add.make_args(E5exp):
    cf, mon = term.as_independent(*allvars, as_Mul=True)
    exp = [0]*18
    for f_ in sp.Mul.make_args(mon):
        if f_.is_Symbol: exp[allvars.index(f_)] += 1
        elif f_.is_Pow and f_.exp.is_integer: exp[allvars.index(f_.base)] += int(f_.exp)
        else: raise ValueError("factor")
    terms.append((tuple(exp), sp.Rational(cf)))
check("B1 E5 parsed: 1910 terms, log-weight degrees <= 5", len(terms) == 1910)
import dod5d_exact as DX
E5X = DX.expand_E5(terms)
stored = pickle.load(open(HERE+"/E5_dod_expanded.pkl","rb"))
same = len(E5X) == len(stored) and all(E5X.get(k) == tuple(v) for k, v in stored.items())
check("B2 exact expansion reproduces E5_dod_expanded.pkl (11006 monomials)", same, len(E5X))
# spot check vs symeng.E5_sym at t=0.37
mp.mp.dps = 50
import symeng
r3 = mp.sqrt(3)
M = [((1 if j<=3 else -1)*mp.cos(j*mp.pi/6), (1 if j<=3 else -1)*mp.sin(j*mp.pi/6)) for j in range(6)]
P = symeng.edge_params(M, (2-r3, mp.mpf(0)), ((2*r3-3)/2, (2-r3)/2))
u = mp.sqrt(2-r3)
KAPm = [mp.mpf(str(k[0]))+mp.mpf(str(k[1]))*r3 for k in DX.KAP]
logCA = mp.log((2-r3)/2); logCC0 = mp.log(2-r3); logCC1 = mp.log((2*r3-3)/2)
def E5_eval(E, t):
    LT = mp.log(t); LU = mp.log(1-t); W = [mp.log(1+KAPm[i]*t) for i in range(4)]
    tot = mp.mpf(0)
    for (n, le), (p, q) in E.items():
        v = t**n * LT**le[0] * LU**le[1]
        for i in range(4):
            if le[2+i]: v *= W[i]**le[2+i]
        if le[6]: v *= logCA**le[6]
        if le[7]: v *= logCC0**le[7]
        if le[8]: v *= logCC1**le[8]
        tot += (mp.mpf(str(p))+mp.mpf(str(q))*r3)*v
    return u*tot
d37 = abs(E5_eval(E5X, mp.mpf('0.37')) - symeng.E5_sym(mp.mpf('0.37'), M, P, conv='B'))
check("B3 exact E5 vs symeng.E5_sym at t=0.37", d37 < mp.mpf('1e-40'), mp.nstr(d37, 3))

# ---- itab spot check + assembly ----
import flint
IT = pickle.load(open(HERE+"/itab140.pkl","rb"))
def one_monomial(m, DPS=60):
    flint.ctx.prec = int(DPS*3.322)+40
    def q3f(c): return flint.arb(str(c[0])) + flint.arb(str(c[1]))*flint.arb(3).sqrt()
    KAPa = [q3f(k) for k in DX.KAP]
    PI = flint.arb.pi(); ONE = flint.arb(1); TWO = flint.arb(2)
    n,a,b,c0,c1,c2,c3 = m
    h = ONE; prev = None
    for level in range(12):
        h = h/2; acc = flint.arb(0); k = 0
        while True:
            u_ = k*h
            arg = (PI/2)*u_.sinh()
            th = arg.tanh()
            tmr = TWO/((TWO*arg).exp()+ONE)
            w = (PI/4)*u_.cosh()/(arg.cosh()**2)*h
            tpl = (ONE+th)/2; om = tmr/2
            ts = [(tpl, om)] if k == 0 else [(tpl, om), (om, tpl)]
            contrib = 0.0
            for tt, omt in ts:
                f_ = tt**n * tt.log()**a * omt.log()**b
                for i, cc in enumerate([c0,c1,c2,c3]):
                    if cc: f_ = f_*(ONE+KAPa[i]*tt).log()**cc
                contrib = max(contrib, abs(float((w*f_).mid())))
                acc = acc + w*f_
            k += 1
            if contrib < 10.0**(-(DPS-8)) and k > 10: break
        if prev is not None and abs(float((acc-prev).mid())) < 10.0**(-(DPS-10)): return acc
        prev = acc
    return acc
random.seed(5)
mons = sorted(IT)
sample = random.sample(mons, 6) + [(0,0,0,0,0,0,0), (0,0,5,0,0,0,0)]
worst = mp.mpf(0)
for m in sample:
    v = one_monomial(m)
    d = abs(mp.mpf(IT[m]) - mp.mpf(v.mid().str(65, radius=False)))
    worst = max(worst, d)
check("B4 itab140 spot check (8 monomials, fresh quadrature)", worst < mp.mpf('1e-50'), mp.nstr(worst, 3))
mp.mp.dps = 150
r3 = mp.sqrt(3); u = mp.sqrt(2-r3)
CA = mp.log((2-r3)/2); CC0 = mp.log(2-r3); CC1 = mp.log((2*r3-3)/2)
P5u = P5_170/u
tot = mp.mpf(0)
for (n, le), coef in E5X.items():
    v = (mp.mpf(str(coef[0]))+mp.mpf(str(coef[1]))*r3)*mp.mpf(IT[(n, le[0], le[1], le[2], le[3], le[4], le[5])])
    if le[6]: v *= CA**le[6]
    if le[7]: v *= CC0**le[7]
    if le[8]: v *= CC1**le[8]
    tot += v
res_itab = abs(1440*tot - P5u)
check("B5 itab assembly of P^(5)(1)/u vs engine (>=60 digits)", res_itab < mp.mpf('1e-100'), mp.nstr(res_itab, 3))

# ---------------- (C) rewrite + G-table ----------------
import dod5d_grewrite as GR
rewrite = pickle.load(open(HERE+"/rewrite_exact.pkl","rb"))
random.seed(7)
samp = random.sample(sorted(rewrite), 60)
ok = all(GR.monomial_integral(*m) == rewrite[m] for m in samp)
check("C1 60 random exact rewrites reproduce rewrite_exact.pkl", ok)
G1raw = pickle.load(open(HERE+"/gtab6_d130_N320.pkl","rb"))
mp.mp.dps = 130
G1 = {k: mp.mpf(v) for k, v in G1raw.items()}
p0 = 4+2*r3; p2 = -(1+r3)/2
c_quad = abs(G1[(2,4)] - mp.quad(lambda t: mp.log(1-t/p2)/(t-p0), [0,1]))
c_shuf = abs(G1[(2,)]*G1[(4,)] - G1[(2,4)] - G1[(4,2)])
c_poly = abs(G1[(0,0,0,0,2)] + mp.polylog(5, 1/p0))
c_log = abs(G1[(5,)] - mp.log(1+(r3-3)/3))
check("C2 G-table self-checks (quad, shuffle, polylog, log)",
      c_quad < mp.mpf('1e-120') and c_shuf < mp.mpf('1e-120') and c_poly < mp.mpf('1e-90') and c_log < mp.mpf('1e-120'),
      [mp.nstr(z,3) for z in (c_quad, c_shuf, c_poly, c_log)])
SLOT = {GR.L0: 0, GR.L1: 1, GR.P0: 2, GR.P1: 3, GR.P2: 4, GR.P3: 5}
mp.mp.dps = 150
tot = mp.mpf(0)
for (n, le), coef in E5X.items():
    rw = rewrite[(n, le[0], le[1], le[2], le[3], le[4], le[5])]
    iv = mp.mpf(0)
    for w, cw in rw.items():
        iv += (mp.mpf(str(cw[0]))+mp.mpf(str(cw[1]))*r3)*G1[tuple(SLOT[x] for x in w)]
    v = (mp.mpf(str(coef[0]))+mp.mpf(str(coef[1]))*r3)*iv
    if le[6]: v *= CA**le[6]
    if le[7]: v *= CC0**le[7]
    if le[8]: v *= CC1**le[8]
    tot += v
res_gtab = abs(1440*tot - P5u)
check("C3 gtab assembly of P^(5)(1)/u vs engine (>=60 digits)", res_gtab < mp.mpf('1e-90'), mp.nstr(res_gtab, 3))

# ---------------- (D) symbols, obstruction, controls ----------------
from qs3_symbol import (q3, q3add, q3mul, Q3V, ZERO, ONE, MONE, add, scale, shuffle,
                        g_symbol, li_symbol, const_pow, logvec)
import obstruction32 as OB
GT = pickle.load(open(HERE+"/gsym_tab.pkl","rb"))
random.seed(9)
wsamp = random.sample(sorted(GT), 24)
ok = all(g_symbol(ZERO, w, ONE) == GT[w] for w in wsamp)
check("D1 24 random g_symbols reproduce gsym_tab.pkl", ok)
# re-assemble total symbol (grouped fast path)
from collections import defaultdict
E5X_ = stored
CA_VEC = (0, -2, 0, 0); CC0_VEC = (-1, 0, 0, 0); CC1_VEC = (0, -2, 1, 0)
A = defaultdict(lambda: {})
for (n, le), coef in E5X_.items():
    if sum(le) != 5: continue
    vwt = le[0]+le[1]+le[2]+le[3]+le[4]+le[5]
    cpat = (le[6], le[7], le[8])
    rw = rewrite[(n, le[0], le[1], le[2], le[3], le[4], le[5])]
    Ag = A[cpat]
    for w, cw in rw.items():
        if len(w) != vwt: continue
        Ag[w] = q3add(Ag.get(w, ZERO), q3mul(coef, cw))
for cpat in A: A[cpat] = {w: c for w, c in A[cpat].items() if c != ZERO}
def const_sym(e6, e7, e8):
    out = {(): ONE}
    if e6: out = shuffle(out, const_pow(CA_VEC, e6))
    if e7: out = shuffle(out, const_pow(CC0_VEC, e7))
    if e8: out = shuffle(out, const_pow(CC1_VEC, e8))
    return out
total = {}
for cpat, Ag in A.items():
    s = {}
    for w, c in Ag.items():
        s = add(s, scale(GT[w], c))
    if any(cpat): s = shuffle(s, const_sym(*cpat))
    total = add(total, s)
total = scale(total, q3(1440))
D = pickle.load(open(HERE+"/dod5d_symbols.pkl","rb"))
check("D2 total symbol re-assembled == dod5d_symbols.pkl (243 entries)", add(total, scale(D['total'], MONE)) == {})
tT, _ = OB.obstruction(D['total'], 3)
check("D3 TOTAL obstruction nonzero (24 pairings, Lie_3(x)Lie_2 Hall bases)",
      not OB.is_zero(tT) and all(v != 0 for row in tT for v in row),
      "e.g. [X,[X,Y]]x[X,Y] = (%s, %s)" % (tT[0][0][0], tT[0][0][1]))
tloc, _ = OB.obstruction(D['blocks']['local'], 3)
check("D4 local block obstruction == 0", OB.is_zero(tloc))
acc = [[0]*3 for _ in range(8)]
for blk in D['blocks']:
    tb, _ = OB.obstruction(D['blocks'][blk], 3)
    for i in range(8):
        for j in range(3): acc[i][j] = acc[i][j] + tb[i][j]
check("D5 sum of block obstructions == total", all(acc[i][j] == tT[i][j] for i in range(8) for j in range(3)))
mir = [('single-0','single-1'), ('single-2','single-3'), ('double-03','double-12'),
       ('double-02','double-13'), ('triple-023','triple-123'), ('triple-012','triple-013')]
okm = True
for a, b in mir:
    ta, _ = OB.obstruction(D['blocks'][a], 3); tb, _ = OB.obstruction(D['blocks'][b], 3)
    okm &= all(ta[i][j] == tb[i][j] for i in range(8) for j in range(3))
check("D6 mirror equalities of block obstructions", okm)
# controls
pts = [q3(2,-1), q3(F(-1,2),F(1,2)), q3(0,F(1,3)), q3(7,-4), q3(1,F(-1,2)),
       q3(-1,F(2,3)), q3(F(1,2),F(-1,6)), q3(F(1,3)), q3(1), q3(F(1,2)),
       q3(-1,1), q3(3,2), q3(5,2)]
okA = all(OB.is_zero(OB.obstruction(li_symbol(5, z), 3)[0]) for z in pts)
check("D7 control: Li_5 at 13 S-units pairs to 0", okA)
random.seed(11); okB = True
for _ in range(300):
    k = random.choice([1, 2])
    w1 = tuple(random.choice([0,1,2]) for _ in range(k))
    w2 = tuple(random.choice([0,1,2]) for _ in range(5-k))
    s1 = {w1: q3(random.randint(-4,4), random.randint(-3,3))}
    s2 = {w2: q3(random.randint(-4,4), random.randint(-3,3))}
    if not OB.is_zero(OB.obstruction(shuffle(s1, s2), 3)[0]): okB = False; break
check("D8 control: 300 random shuffle products pair to 0", okB)
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2 as Q
symj = json.loads(open(THEORY+"M4c_w5_closed_form/symbol.json").read())['symbol']
Sm4 = {tuple(map(int, w)): Q(F(p), F(q)) for w, (p, q) in symj.items()}
tm4, _ = OB.obstruction(Sm4, 2)
check("D9 control: M4c W5 pairs to (84, -38 sqrt2)", tm4[0][0] == Q(84) and tm4[1][0] == Q(0,-38))
bs = pickle.load(open(THEORY+"OCT5D2_octagon_depth_test/block_symbols.pkl","rb"))
def loadq2(d): return {tuple(map(int, w)): Q(F(p), F(q)) for w, (p, q) in d.items()}
toct, _ = OB.obstruction(loadq2(bs['total']), 2)
tloc8, _ = OB.obstruction(loadq2(bs['blocks']['local']), 2)
check("D10 control: octagon total (-1485+1035 sqrt2, 900-540 sqrt2), local 0",
      toct[0][0] == Q(-1485,1035) and toct[1][0] == Q(900,-540) and OB.is_zero(tloc8))
# G-symbol identities
s1 = g_symbol(ZERO, (ZERO, q3(4,2)), ONE)
from qs3_symbol import q3div
s3 = g_symbol(ZERO, (ZERO,)*4+(q3(4,2),), ONE)
ok1 = add(s1, li_symbol(2, q3div(ONE, q3(4,2)))) == {}
ok3 = add(s3, li_symbol(5, q3div(ONE, q3(4,2)))) == {}
s2a = shuffle(g_symbol(ZERO, (q3(4,2),), ONE), g_symbol(ZERO, (q3(-3,-2),), ONE))
s2b = add(g_symbol(ZERO, (q3(4,2), q3(-3,-2)), ONE), g_symbol(ZERO, (q3(-3,-2), q3(4,2)), ONE))
ok2 = add(s2a, scale(s2b, MONE)) == {}
check("D11 G-symbol identities (0,p0)=-Li_2(1/p0), (0^4,p0)=-Li_5(1/p0), shuffle(p0,p1)", ok1 and ok2 and ok3)
# Li_3,2 correction
sol = pickle.load(open(HERE+"/li32_solution.pkl","rb"))
L3b = OB.lie3_basis(3); L2b = OB.lie2_basis(3)
combo = [Q3V(0)]*24
for ab, xx, yy, c in sol:
    w = (ZERO,)*(ab[1]-1) + (q3div(ONE, yy),) + (ZERO,)*(ab[0]-1) + (q3div(ONE, q3mul(xx, yy)),)
    S = g_symbol(ZERO, w, ONE)
    Psi = OB.psi_tensor(S)
    for ii, (_, AA) in enumerate(L3b):
        for jj, (_, BB) in enumerate(L2b):
            combo[ii*3+jj] = combo[ii*3+jj] + c*Q3V(OB.pair(Psi, AA, BB))
resid = [combo[ii*3+jj] - Q3V(tT[ii][jj]) for ii in range(8) for jj in range(3)]
check("D12 24-term Li_3,2/Li_4,1 correction == total obstruction (residual exactly 0)",
      all(r == 0 for r in resid))

print("\\ntotal time %.1fs" % (time.time()-T0))
if failures:
    print("FAILURES:", failures); sys.exit(1)
print("ALL CHECKS PASS — depth two needed (PROVED symbol level).")
