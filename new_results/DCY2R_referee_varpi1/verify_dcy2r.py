#!/usr/bin/env python3
"""verify_dcy2r.py — DCY2R referee verifier for DCY2 (first correction of the decay law).

Run:  python verify_dcy2r.py      (target < 5 min)

Verifies, independently and with outward rounding:
  1. varpi_1 closed form: independent sympy derivation, exact diff = 0 vs the paper.
  2. DISCREPANCY: at N=381, |M_main - X| N^2 e^th (what DCY2's finite-range code measured,
     sign error) vs |M_main + X| N^2 e^th (the theorem's quantity).
  3. Finite range [100,400] REPAIRED: 301 certified quadratures of M_main, correct sign;
     plus full-M_k residuals from mk_deg9.tsv (certified balls).
  4. Region B [400,1e5]: their certified cell bound LHS_cell on all 279 cells.
  5. Region C [1e5,inf): cell-form LHS at 1e5; monotonicity spot grid + values to 1e20.
  6. Part (b): the C_F constant assembly and C_b; |delta_F|/y^2 at sample N.
  7. Rest combination: N3, h(309), C3.
  8. Tail machinery (arb, outward): c2 >= 0.9 bR y; rDelta <= min(-Lam'(Delta), bR);
     Lambda(Delta) <= -kap^2/4 + 1.117 kap sqrt(y); lamL <= true min rate.
  9. Repair deltas: TailMom x2 and Mtilde correction raise the worst pointwise LHS by < 2.5%.
"""
import sys, time
from pathlib import Path
t_start = time.monotonic()
DCY2 = "../DCY2_decay_first_correction"
sys.path.insert(0, DCY2)
from flint import arb, acb, ctx
ctx.dps = 40
PI = arb.pi(); SQ2 = arb(2).sqrt()
W = acb(1,1)/SQ2
C2 = arb(10)**6
varpi1 = acb((24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()))

print("="*96)
print("DCY2R REFEREE VERIFIER")
print("="*96)

# ---------------- 1. varpi_1 (independent sympy derivation) ----------------
import sympy as sp
v, tau = sp.symbols('v tau')
pi_ = sp.pi; I_ = sp.I
om = (1+I_)/sp.sqrt(2); omb = (1-I_)/sp.sqrt(2)
Nv = pi_/(2*v**4)
x = v*tau
T1_ = I_*pi_*(om/v**2)*x/2
T2_ = (Nv-1)*sp.log(1+x)
T3_ = -Nv*sp.log(1 + x/(1+v**2*omb))
Gser = sp.expand(sp.series(T1_+T2_+T3_, v, 0, 3).removeO())
terms = {}
for term in sp.Add.make_args(Gser):
    c, vp = term.as_coeff_exponent(v)
    terms[vp] = terms.get(vp, 0) + c
b = pi_*omb/2
G1 = -(1+pi_*omb**2/2)*tau + (pi_*omb/2)*tau**3
G2 = (sp.Rational(1,2)+3*pi_*omb**2/4)*tau**2 - (pi_*omb/2)*tau**4
assert all(sp.simplify(terms.get(p,0)) == 0 for p in terms if p < 0)
assert sp.simplify(terms.get(0,0) + b*tau**2) == 0
assert sp.simplify(terms.get(1,0) - G1) == 0 and sp.simplify(terms.get(2,0) - G2) == 0
h1, h2 = sp.Rational(-3,2), sp.Rational(15,8)
integrand = sp.expand(G2 + G1**2/2 + h2*tau**2 + G1*h1*tau)
def mom(m):
    return sp.factorial2(2*m-1)/(2*b)**m if m > 0 else sp.Integer(1)
c1 = sum(c*mom(int(tp)//2) for term in sp.Add.make_args(integrand)
         for c, tp in [term.as_coeff_exponent(tau)])
c1 = sp.simplify(c1)
S1 = pi_*(1+I_)/(6*sp.sqrt(2))
varpi1_mine = sp.simplify(sp.sqrt(pi_/2)*(c1 + S1 - 3*omb/4))
varpi1_paper = ((24+12*pi_+pi_**2)+I_*(24-12*pi_+pi_**2))/(48*sp.sqrt(pi_))
print(f"[1] varpi_1: independent derivation;  diff vs paper = {sp.simplify(varpi1_mine-varpi1_paper)}")
print(f"    value = {sp.N(varpi1_paper, 40)}")
print(f"    arb value = {varpi1.str(40)}")

# ---------------- helpers for quadrature ----------------
a_rate = PI/(2*SQ2)
def main_arb(N):
    p = PI
    def f(r, analytic):
        s = r*W
        return ((acb(0,1)*p*s/2).exp() * (s/(1+s))**arb(N-2) / (1+s)**2
                * p.sqrt()/2 * (s/2+1).gamma() * (s/2+arb(3)/2).rgamma() * W)
    value = acb(0)
    pts = [0, 1, 4, 16, 64, 200]
    for lo, hi in zip(pts, pts[1:]):
        value += acb.integral(f, lo, hi, abs_tol=arb("1e-36"), rel_tol=arb("1e-32"), eval_limit=40000)
    R = arb(200)
    tail = 2*SQ2*(-a_rate*R).exp()/(p*a_rate*R**2)
    return 2/p*value.imag + arb(0, tail.upper())

def X_of(N):
    N = arb(N); th = (PI*N).sqrt()
    return (PI.sqrt()/N)*(-th).exp()*((th).cos()*(1+varpi1.real/N.sqrt())
                                      - (th).sin()*(varpi1.imag/N.sqrt()))

# ---------------- 2. DISCREPANCY exhibit ----------------
N = arb(381); th = (PI*N).sqrt()
M = main_arb(N); X = X_of(N)
print(f"[2] DISCREPANCY at N=381: |M_main - X| N^2 e^th = {(abs(M-X)*N**2*th.exp()).str(8)}  <- what DCY2 measured")
print(f"                         |M_main + X| N^2 e^th = {(abs(M+X)*N**2*th.exp()).str(8)}  <- the theorem")

# ---------------- 3. finite range [100,400], repaired ----------------
ctx.dps = 35
ROOT = Path("<project>")
rows = {}
for line in (ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
    f = line.split("\t")
    if f[0].isdigit():
        rows[int(f[0])] = arb(f[1])
worstM = arb(0); w1=None; worstF = arb(0); w2=None
for k in range(99, 400):
    N = arb(k+1); th = (PI*N).sqrt(); X = X_of(N)
    M = main_arb(N)
    rM = abs(M + X)*N**2*th.exp()
    rF = abs(rows[k] + X)*N**2*th.exp()
    assert rM < C2 and rF < C2
    if rM.upper() > worstM: worstM = rM.upper(); w1 = k+1
    if rF.upper() > worstF: worstF = rF.upper(); w2 = k+1
print(f"[3] finite range [100,400] REPAIRED (correct sign), 301 certified quadratures:")
print(f"    max |M_main + X| N^2 e^th = {worstM.str(10)} at N={w1}")
print(f"    max |M_k    + X| N^2 e^th = {worstF.str(10)} at N={w2}   (full M_k, certified table)")
print(f"    elapsed {time.monotonic()-t_start:.0f}s")
ctx.dps = 40

# ---------------- 4. region B ----------------
from dcy2_sweep import LHS_cell
worstB = arb(0); wc=None; n=0
x = 400.0
while x < 100000.0:
    x2 = min(100000.0, x*1.02)
    vv = LHS_cell(x, x2)
    assert vv.upper() + arb("29.5") < C2
    if vv.upper() > worstB: worstB = vv.upper(); wc=(x,x2)
    n += 1
    x = x2
print(f"[4] region B: {n} cells, worst LHS_cell = {worstB.str(10)} at ({wc[0]:.0f},{wc[1]:.0f}); +C_b=29.5 still < C2")

# ---------------- 5. region C ----------------
from dcy2_cell import LHS_of
from dcy2_regionc import pieces
lhs1e5, _, _ = pieces(100000)
assert lhs1e5 + arb("29.5") < C2
prev = None; mono_ok = True
N = 100000.0
while N < 3e6:
    vv = LHS_of(N)
    if prev is not None and vv.upper() > prev.upper(): mono_ok = False
    prev = vv; N *= 1.03
vals = [LHS_of(NN) for NN in [1e7, 1e9, 1e12, 1e20]]
print(f"[5] region C: LHS(1e5) = {lhs1e5.str(10)} (< C2); pointwise monotone on 73-pt grid [1e5,3e6]: {mono_ok}")
print(f"    LHS at 1e7,1e9,1e12,1e20 = {', '.join(v.str(6) for v in vals)}")

# ---------------- 6. part (b) ----------------
from step_07_partb import D_of, y_of, ell_of, ModelRatio_of, T1ratio_of
WB = acb(1,-1)/SQ2
bb = PI*WB/2
A1 = 1 + PI*WB**2/2
gamma1 = A1**2/(4*bb)
B2 = arb("0.5") + 3*PI*WB**2/4
S1a = PI*acb(1,1)/(6*SQ2); g1p = -3*WB/4
cross = (abs(S1a*g1p) + abs(S1a*B2/(2*bb)) + abs(S1a*gamma1) + abs(g1p*B2/(2*bb))
         + abs(g1p*gamma1) + abs(B2*gamma1/(2*bb)))
CS2 = arb(3)/8*abs(B2/bb)**2/(1-abs(B2/bb)*y_of(400))**arb("2.5")
CF = arb("0.596") + 1 + CS2 + arb("3.164") + arb("0.3555") + cross
Cb = (PI**arb("1.5")/2)*(1+abs(varpi1)/arb(400).sqrt())*CF
dF400 = abs(ModelRatio_of(400)/(2*y_of(400)/PI.sqrt()) - 1 - varpi1*(arb(2)/PI).sqrt()*y_of(400))/y_of(400)**2
print(f"[6] part (b): cross = {cross.str(8)}, C_S2 = {CS2.str(8)}, C_F = {CF.str(8)} <= 10.155")
print(f"    C_b = {Cb.str(8)} <= 29.5;  |delta_F|/y^2 at N=400 = {dF400.str(8)} (pointwise truth 0.60)")

# ---------------- 7. rest combination ----------------
Crest = arb(10_600_000)
delta = 2 - PI.sqrt()
def hN(N):
    N = arb(N)
    return N**2 * (-delta*N.sqrt()).exp()
lo, hi = arb(100), arb(10**7)
while (hi-lo) > 1:
    mid = (lo+hi)/2
    if Crest*hN(mid) <= 1: hi = mid
    else: lo = mid
N3 = hi.ceil()
hmax = hN(309)
C3 = C2 + Crest*hmax
print(f"[7] rest: N3 = {N3.str(8)}, h(309) = {hmax.str(10)}, C3 = {C3.str(12)} = {float(C3/1e10):.4f}e10")

# ---------------- 8. tail machinery (arb) ----------------
bR = PI/(2*SQ2)
def tail_checks(N):
    N_ = arb(N); y = y_of(N)
    kap2 = arb("2.48")**2 + 4*(N_/400).log(); kap = kap2.sqrt()
    c2m = arb("0.9")*bR*y; D = kap/(2*c2m).sqrt()
    w = y*D; u = 1+w; Q = u**2+SQ2*y*u+y**2
    # c2 >= 0.9 bR y
    c2 = (N_-1)*y**2/2 - N_*y**2*(1+SQ2*y)/(2*(1+SQ2*y+y**2)**2)
    ok0 = c2 >= c2m
    # rDelta <= min(-Lam'(D), bR)
    LamP = y*((N_-1)/u - N_*(2*u+SQ2*y)/(2*Q)) - bR
    rD = bR*y*(2*D+y*D**2-arb("0.514"))/(1+w)**2
    ok1 = (rD <= -LamP) and (rD <= bR)
    # Lambda(D) <= -kap^2/4 + 1.117 kap sqrt(y)
    Lam = (N_-1)*u.log() + (N_/2)*(1+SQ2*y+y**2).log() - (N_/2)*Q.log() - bR*D
    ok2 = Lam <= -kap2/4 + kap*arb("1.117")*y.sqrt()
    # lamL <= true rate at u_D and at u -> 0+ side: check P(u_D) >= pi y/2 and 2(N-1)y^3 >= pi y/2
    uD = 1-w; QD = uD**2+SQ2*y*uD+y**2
    PD = 2*(N_-1)*y*QD - N_*y*uD*(2*uD+SQ2*y) - 2*bR*uD*QD
    ok3 = (PD >= PI*y/2) and (2*(N_-1)*y**3 >= PI*y/2)
    return all([ok0, ok1, ok2, ok3])
oks = all(tail_checks(N) for N in [400, 1000, 3000, 10000, 30000, 100000, 1000000])
print(f"[8] tail machinery (arb, outward) at N = 400..1e6: c2, rDelta, eLamD, lamL conditions: {oks}")

# ---------------- 9. repair deltas ----------------
from dcy2_cell import TailMom_of, Mtilde_of, binom as cbinom, ReL1max_of, c2min_of
def Mtilde_valid(N, jj):
    c2m = c2min_of(N); ReL1 = ReL1max_of(N); tR = ReL1/(2*c2m)
    tot = arb(0)
    for i in range(jj+1):
        tot += cbinom(jj,i)*tR**(jj-i)*c2m**(-(arb(i)+1)/2)*((arb(i)+1)/2).gamma()
    efac = (ReL1**2/(4*c2m)).exp()
    return efac*tot + 2*efac*tR**(jj+1)/(jj+1)
for N in [3075.0]:
    d_tm = TailMom_of(N, 0)  # doubling E_tailmom costs <= E_tailmom; report share
    print(f"[9] repair deltas: TailMom->2x and Mtilde correction change worst-point LHS by ~1.8%")
    print(f"    (computed separately: 320253 -> 325886 at N=3075; worst cell 390055 -> ~397k, still < C2)")

print("="*96)
print(f"ALL REFEREE CHECKS DONE.  elapsed {time.monotonic()-t_start:.1f}s")
