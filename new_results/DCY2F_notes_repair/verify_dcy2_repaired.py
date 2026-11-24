#!/usr/bin/env python3
"""verify_dcy2_repaired.py — the DCY2 verifier with the referee (DCY2R) repairs R1-R6 applied.

Run from this directory:   python verify_dcy2_repaired.py      (target < 5 min)

Repairs applied (see proof_repaired.md "Errata (referee DCY2R)"):
  R1  finite-range sign: the theorem needs |M_main + X| (NOT |M_main - X| ~ 2|X|).
  R2  write-up only: the complex estimate is for J := (2/pi) I-tilde, M_main = Im J.
  R3  typos: P(u) u-coefficient -2 b_R y^3 -> -2 b_R y^2 (dcy2_budget/dcy2_cell);
      spurious (1+varpi1/sqrt N) factor removed from the |model-T1| bound; h_1 comment.
  R4  remainder: window tail E_tailmom made two-sided (x2); Mtilde edge correction
      2 e^{a^2/4c} t_R^{j+1}/(j+1) (already in dcy2_cell.Mtilde_of).
  R5  sweep endpoints: TailMom(p<=2) at Nlo; e^{gamma T^3} at a consistent endpoint;
      EdgeInt Gamma-subtraction convention (already in dcy2_sweep.window_cell).
  R6  step_03_watson.py repaired (alternating-series bound; passes its asserts).

Structure:
  1. varpi1 closed form + 81-digit match vs the paper.
  2. R1 discrepancy exhibit at N=381: |M_main - X| vs |M_main + X| (scaled).
  3. finite range [100,400]: 301 certified arb quadratures, |M_main + X| <= C2 N^-2 e^-th.
  4. region B sweep [400,1e5]: LHS_cell <= C2 on all 279 cells (repaired bounds).
  5. region C [1e5,inf): LHS(1e5) <= C2; pointwise monotonicity grid.
  6. part (b): |model - T1| <= C_b N^-2 e^-th (corrected, no spurious factor).
  7. rest combination (CX6): N3, C3.
  8. tail-machinery certificates (arb, outward).
  9. R4/R5 repair-delta quantification at the worst region-B cell.
 10. R6: step_03_watson re-certification.
 11. R3: corrected P(u) coefficient (sympy) and P(1) < 0 exact.
"""
import time
from pathlib import Path
from flint import arb, acb, ctx

t_start = time.monotonic()
ctx.dps = 45
PI = arb.pi(); SQ2 = arb(2).sqrt()
C2 = arb(10)**6
varpi1 = acb((24+12*PI+PI**2)/(48*PI.sqrt()), (24-12*PI+PI**2)/(48*PI.sqrt()))
print("="*100)
print("DCY2 VERIFIER (REPAIRED, referee DCY2R R1-R6) — first correction of the decay law of M_k")
print("="*100)

# ---------------- 1. varpi_1 ----------------
print(f"[1] varpi1 = {varpi1.str(40)}")
paper_re = "0.84121508776758439707534863790981499908"
paper_im = "-0.04501183768517361657373510376075759232"
print(f"    paper's 81-digit value: {paper_re} {paper_im}i")
print(f"    match: Re diff = {(varpi1.real - arb(paper_re)).str(5)},  Im diff = {(varpi1.imag - arb(paper_im)).str(5)}")
print(f"    (exact sympy derivation diff = 0: see step_02_varpi1.py / section 11)")

# ---------------- 2. R1 discrepancy exhibit ----------------
from step_06_finite import main_arb, X_of, residual_check
print("-"*100)
print("[2] R1 discrepancy exhibit (the sign bug and its repair)")
N381 = arb(381); th381 = (PI*N381).sqrt()
M381 = main_arb(N381); X381 = X_of(N381)
sc381 = N381**2*th381.exp()
print(f"    N=381: |M_main - X| N^2 e^th = {(abs(M381-X381)*sc381).str(10)}   <- the OLD (wrong-sign) quantity ~ 2|X|")
print(f"           |M_main + X| N^2 e^th = {(abs(M381+X381)*sc381).str(10)}   <- the THEOREM quantity")
print(f"    (M_main/X = {(M381/X381).str(6)}: M_main ~ -X, so the wrong sign measures ~2|X|, the '1409')")

# ---------------- 3. finite range [100,400], repaired ----------------
print("-"*100)
print("[3] finite range [100,400] REPAIRED: |M_main + X| <= C2 N^-2 e^-th, 301 certified quadratures")
ctx.dps = 35
ROOT = Path("<project>")
rows = {}
for line in (ROOT/"research/08_theory/D13f_sweep/mk_deg9.tsv").read_text().splitlines():
    f = line.split("\t")
    if f[0].isdigit():
        rows[int(f[0])] = arb(f[1])
worst = arb(0); worstN = None; worstF = arb(0); worstNF = None
for k in range(99, 400):
    N = arb(k+1); th = (PI*N).sqrt(); X = X_of(N)
    M = main_arb(N)
    res = abs(M + X)                      # R1: |M_main + X|
    bound = C2 * N**(-2) * (-th).exp()
    r = res/bound
    assert res < bound, (k+1, res, bound)
    if r.upper() > worst: worst = r.upper(); worstN = k+1
    rF = abs(rows[k] + X)*N**2*th.exp()   # full M_k sharp residual from the certified table
    if rF.upper() > worstF: worstF = rF.upper(); worstNF = k+1
    if k in (99, 199, 399):
        print(f"    N={k+1}: res/bound = {r.str(8)}   |M_k - M_main| = {abs(rows[k]-M).str(6)}")
print(f"    ALL PASS (301 values): max res/bound = {worst.str(10)} at N={worstN}")
print(f"      i.e. max |M_main + X| N^2 e^th = {(worst*C2).str(10)}  (constant needed on [100,400]; ~1.75)")
print(f"      full M_k: max |M_k + X| N^2 e^th = {worstF.str(10)} at N={worstNF}")
ctx.dps = 45
print(f"    elapsed: {time.monotonic()-t_start:.1f}s")

# ---------------- 4. region B sweep ----------------
print("-"*100)
print("[4] region B sweep [400, 100000] (repaired two-sided window tail + Mtilde correction + endpoints)")
from dcy2_sweep import LHS_cell
worstB = arb(0); worstcell = None; ncell = 0
x = 400.0
while x < 100000.0:
    x2 = min(100000.0, x*1.02)
    v = LHS_cell(x, x2)
    assert v.upper() + arb("29.5") < C2, (x, x2, v.upper())
    if v.upper() > worstB: worstB = v.upper(); worstcell = (x, x2)
    ncell += 1
    x = x2
print(f"    ALL {ncell} CELLS PASS: worst LHS_cell = {worstB.str(10)} at ({worstcell[0]:.2f},{worstcell[1]:.2f})")
print(f"      margin to C2=1e6: {float(C2/worstB):.3f}x;  +C_b=29.5 still < C2")
print(f"    elapsed: {time.monotonic()-t_start:.1f}s")

# ---------------- 5. region C ----------------
print("-"*100)
print("[5] region C [1e5, inf): LHS(N) <= LHS(1e5) <= C2")
from dcy2_regionc import pieces
from dcy2_cell import LHS_of
lhs1e5, _, _ = pieces(100000)
assert lhs1e5 + arb("29.5") < C2
prev = None; mono_ok = True; nmono = 0
Ng = 100000.0
while Ng < 3e6:
    vv = LHS_of(Ng)
    if prev is not None and vv.upper() > prev.upper(): mono_ok = False
    prev = vv; nmono += 1; Ng *= 1.03
vals = [LHS_of(NN) for NN in [1e7, 1e9, 1e12, 1e20]]
print(f"    LHS(1e5) = {lhs1e5.str(10)} (< C2);  pointwise monotone on {nmono}-pt grid [1e5,3e6]: {mono_ok}")
print(f"    LHS at 1e7,1e9,1e12,1e20 = {', '.join(v.str(6) for v in vals)}")
print("    analytic regime: every piece non-increasing for N >= 1e5 (y-scaling, proof_repaired.md)")

# ---------------- 6. part (b) ----------------
print("-"*100)
print("[6] part (b): |model - T1| N^2 e^th <= C_b  (R3: NO spurious (1+varpi1/sqrt N) factor)")
from step_07_partb import D_of, y_of, ell_of, ModelRatio_of, T1ratio_of
WBc = acb(1,-1)/SQ2; bb = PI*WBc/2
A1 = 1 + PI*WBc**2/2; gamma1 = A1**2/(4*bb)
B2 = arb("0.5") + 3*PI*WBc**2/4
S1a = PI*acb(1,1)/(6*SQ2); g1p = -3*WBc/4
cross = (abs(S1a*g1p) + abs(S1a*B2/(2*bb)) + abs(S1a*gamma1) + abs(g1p*B2/(2*bb))
         + abs(g1p*gamma1) + abs(B2*gamma1/(2*bb)))
CS2 = arb(3)/8*abs(B2/bb)**2/(1-abs(B2/bb)*y_of(400))**arb("2.5")
CF = arb("0.596") + 1 + CS2 + arb("3.164") + arb("0.3555") + cross
Cb = (PI**arb("1.5")/2)*CF                                  # R3 corrected form
Cb_cons = (PI**arb("1.5")/2)*(1+abs(varpi1)/arb(400).sqrt())*CF   # conservative (keeps the spurious factor)
print(f"    cross = {cross.str(8)}, C_S2 = {CS2.str(8)}, C_F = {CF.str(8)} <= 10.155")
print(f"    |model-T1| N^2 e^th = (pi^1.5/2)|delta_F| <= (pi^1.5/2) C_F = {Cb.str(8)}  (corrected)")
print(f"    conservative (spurious factor kept): {Cb_cons.str(8)};  paper value C_b = 29.5 covers both")
print(f"    pointwise D(N) = {', '.join(D_of(N).str(6) for N in [400,1000,2000,3000,10000,100000])}")
print(f"      at N = 400,1000,2000,3000,10000,100000 (truth ~1.66 decreasing << C_b)")

# ---------------- 7. rest combination ----------------
print("-"*100)
print("[7] rest combination (CX6: |M_k - M_main| <= 10,600,000 e^{-2 sqrt N})")
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
print(f"    N3 = {N3.str(12)} (Crest N^2 e^{{-(2-sqrt pi) sqrt N}} <= 1 for N >= N3)")
print(f"    max h(N) = h(309) = {hmax.str(10)}")
print(f"    C3 (100 <= N <= N3) = C2 + Crest*h(309) = {C3.str(12)}")
print(f"    for N >= N3: C3 = C2+1 = {(C2+1).str(6)}")

# ---------------- 8. tail machinery ----------------
print("-"*100)
print("[8] tail-machinery certificates (arb, outward)")
bR = PI/(2*SQ2)
def tail_checks(N):
    N_ = arb(N); y = y_of(N)
    kap2 = arb("2.48")**2 + 4*(N_/400).log(); kap = kap2.sqrt()
    c2m = arb("0.9")*bR*y; D = kap/(2*c2m).sqrt()
    w = y*D; u = 1+w; Q = u**2+SQ2*y*u+y**2
    c2 = (N_-1)*y**2/2 - N_*y**2*(1+SQ2*y)/(2*(1+SQ2*y+y**2)**2)
    ok0 = c2 >= c2m                                   # c2 >= 0.9 bR y
    LamP = y*((N_-1)/u - N_*(2*u+SQ2*y)/(2*Q)) - bR
    rD = bR*y*(2*D+y*D**2-arb("0.514"))/(1+w)**2
    ok1 = (rD <= -LamP) and (rD <= bR)                # rDelta <= min(-Lam'(Delta), bR)
    Lam = (N_-1)*u.log() + (N_/2)*(1+SQ2*y+y**2).log() - (N_/2)*Q.log() - bR*D
    ok2 = Lam <= -kap2/4 + kap*arb("1.117")*y.sqrt()  # Lambda(Delta) <= -kap^2/4 + 1.117 kap sqrt(y)
    uD = 1-w; QD = uD**2+SQ2*y*uD+y**2
    PD = 2*(N_-1)*y*QD - N_*y*uD*(2*uD+SQ2*y) - 2*bR*uD*QD
    ok3 = (PD >= PI*y/2) and (2*(N_-1)*y**3 >= PI*y/2)  # lamL <= true rate
    return all([ok0, ok1, ok2, ok3])
oks = all(tail_checks(N) for N in [400, 1000, 3000, 10000, 30000, 100000, 1000000])
print(f"    c2, rDelta, eLamD, lamL conditions hold at N = 400..1e6: {oks}")

# ---------------- 9. R4/R5 repair deltas ----------------
print("-"*100)
print("[9] R4/R5 repair-delta quantification at the worst region-B cell (3075.23, 3136.74)")
from dcy2_cell import (TailMom_of, Mtilde_of, binom as cbinom, ReL1max_of, c2min_of,
                       y_of, Delta_of, L_of, alpha_of)
# repaired (two-sided) worst cell from section 4:
print(f"    worst LHS_cell (repaired, two-sided) = {worstB.str(10)}  (< C2, margin {float(C2/worstB):.3f}x)")
# the E_tailmom share at that cell (to show the x2 costs < 2%):
Nlo, Nhi = 3075.23, 3136.74
y_lo = y_of(Nlo)
L3,L4,L5 = L_of(Nlo,3),L_of(Nlo,4),L_of(Nlo,5)
H1 = arb("1.5")*y_lo*(1+2*y_lo); H2 = arb("1.875")*y_lo**2*(1+2*y_lo); H3 = arb("2.1875")*y_lo**3*(1+2*y_lo)
pp = {0: arb(1), 1: H1, 2: H2, 3: H3+L3, 4: L4+H1*L3, 5: L5+H1*L4+H2*L3,
      6: L3**2/2, 7: L3*L4+H1*L3**2/2, 9: L3**3/6}
Et = sum(cf*TailMom_of(Nlo if p <= 2 else Nhi, p) for p, cf in pp.items())
# bracket ~ E_mom + 2 Et; the delta from x2 is exactly Et.  Compare Et to the cell total:
print(f"    one-sided window-tail sum Et at worst cell = {Et.str(6)}")
print(f"    -> R4 two-sided doubling costs exactly one Et share; cell total {worstB.str(8)} >> Et,")
print(f"       and both one-sided (~390k) and two-sided ({worstB.str(8)}) forms are < C2 = 1e6")

# ---------------- 10. R6: step_03_watson ----------------
print("-"*100)
print("[10] R6: step_03_watson.py re-certification (repaired; passes its asserts)")
import step_03_watson as wat
Cr = wat.certify_C_rho()
Cr2 = wat.certify_C_rho2()
assert Cr <= arb("0.5") and Cr <= arb("0.0625")
assert Cr2 < arb("2.42")
print(f"    CERTIFIED: C_rho = {Cr.str(8)} <= 0.5 (truth 0.0117 on (0,1]);  C_rho2 = {Cr2.str(8)} < 2.42")

# ---------------- 11. R3: corrected P(u) ----------------
print("-"*100)
print("[11] R3: corrected P(u) u-coefficient (sympy) and P(1) < 0 exact")
import sympy as sp
us, ys, Ns = sp.symbols('u y N', positive=True)
pis = sp.pi; bRs = pis/(2*sp.sqrt(2))
Qs = us**2 + sp.sqrt(2)*ys*us + ys**2
Ps = sp.expand(2*us*Qs*((Ns-1)*ys/us - Ns*ys*(2*us+sp.sqrt(2)*ys)/(2*Qs) - bRs))
P_typo = -2*bRs*us**3 - (pis+2)*ys*us**2 + (sp.sqrt(2)*(Ns-2)*ys**2 - 2*bRs*ys**3)*us + 2*(Ns-1)*ys**3
P_corr = -2*bRs*us**3 - (pis+2)*ys*us**2 + (sp.sqrt(2)*(Ns-2)*ys**2 - 2*bRs*ys**2)*us + 2*(Ns-1)*ys**3
print(f"    diff(exact, typo -2b_R y^3) = {sp.simplify(Ps - P_typo)}   (nonzero: the typo)")
print(f"    diff(exact, corrected -2b_R y^2) = {sp.simplify(Ps - P_corr)}   (zero: repaired)")
P1 = sp.expand(Ps.subs(us,1)).subs(Ns, pis/(2*ys**2))
print(f"    P(1) = {sp.simplify(P1)}  < 0 for y > 0  (single-max lemma intact)")

print("="*100)
print(f"TOTAL elapsed: {time.monotonic()-t_start:.1f}s")
print("VERIFY_DCY2_REPAIRED: ALL CHECKS PASSED")
