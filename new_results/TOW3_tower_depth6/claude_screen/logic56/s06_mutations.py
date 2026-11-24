# s06: mutation tests -- perturb one coefficient / constant / sign; the pointwise T1 test MUST fail (except the E3-equivalence control)
import sys, time
sys.dont_write_bytecode = True
from candlib import *
t0 = time.time()
d, m = T3.d, T3.m
I = sp.I
Pi = T3.Pi
base = base_exprs()
DL = [5, 6, 7, 8, 9]
def classify(bad):
    out = {}
    for (p, q, k, delta) in bad:
        typ = 'diag' if p == q else ('adj' if q == p-1 else 'nonadj')
        out[(typ, delta)] = out.get((typ, delta), 0)+1
    return dict(sorted(out.items()))
def run(name, exprs=None, cand_kw=None, t1_kw=None, expect_fail=True, extra=None):
    ex = dict(base)
    if exprs: ex.update(exprs)
    c = Cand(ex, **(cand_kw or {}))
    tested, bad, per = t1_test(c, DL, **(t1_kw or {}))
    status = ('DETECTED' if bad else 'NOT DETECTED') if expect_fail else ('PASS(no failure)' if not bad else 'UNEXPECTED FAILURE')
    print(f'[{status}] {name}: cells {tested}, bad {len(bad)}  by (type,delta): {classify(bad)}  ({round(time.time()-t0,1)}s)', flush=True)
    return len(bad)
# control
run('CONTROL unmutated', expect_fail=False)
# fitted coefficients
run('M1 e(5,1): coeff of d*m  12 -> 13', {(5,1): base[(5,1)] + d*m/(12*Pi(2))})
run('M2 e(5,1): constant -552i -> -552i+1', {(5,1): base[(5,1)] + 1/(12*Pi(2))})
run('M3 e(5,1): coeff of d^4  -4i -> -4i+1', {(5,1): base[(5,1)] + d**4/(12*Pi(2))})
run('M4 e(6,1): coeff of m^2  48 -> 49', {(6,1): base[(6,1)] + m**2/(96*Pi(3))})
run('M5 e(6,1): constant -43398 -> -43397', {(6,1): base[(6,1)] + 1/(96*Pi(3))})
run('M6 e(6,1): coeff of d^6  -16 -> -15', {(6,1): base[(6,1)] + d**6/(96*Pi(3))})
run('M7 e(6,2): constant 268 -> 269', {(6,2): base[(6,2)] + 1/(2*Pi(3))})
run('M8 e(6,2): -4 i m -> -5 i m', {(6,2): base[(6,2)] - I*m/(2*Pi(3))})
run('M9 e(6,2): coeff of d^3  -4 -> -3', {(6,2): base[(6,2)] + d**3/(2*Pi(3))})
# zigzag tops
run('M10 e(5,2) top: A_3=2 -> 3', {(5,2): base[(5,2)]*sp.Rational(3,2)})
run('M11 e(6,3) top: A_4=5 -> 4', {(6,3): base[(6,3)]*sp.Rational(4,5)})
run('M12 e(6,3) top: A_4=5 -> 6', {(6,3): base[(6,3)]*sp.Rational(6,5)})
run('M13 e(5,2) top: factor (d-4) -> (d-3)', {(5,2): base[(5,2)]*(d-3)/(d-4)})
# rho closed forms
run('M14 a_5: numerator constant 25200i -> 25200i+1', {(5,0): base[(5,0)] + 1/(480*(2*d-5)*(2*d-3))})
run('M15 a_6: numerator constant 1984500 -> 1984501', {(6,0): base[(6,0)] + 1/(5760*(2*d-7)*(2*d-5)*(2*d-3))})
run('M16 a_5: coeff of d m^2 (240i) -> 241i', {(5,0): base[(5,0)] + I*d*m**2/(480*(2*d-5)*(2*d-3))})
run('M17 a_6: coeff of m^3 (960i) -> 961i', {(6,0): base[(6,0)] + I*m**3/(5760*(2*d-7)*(2*d-5)*(2*d-3))})
# lower depth sanity (the framework must detect errors at depth 4 and 3 as well)
run('M18 e(4,1): constant -29 -> -28 (TOW2 form)', {(4,1): base[(4,1)] + 1/Pi(2)})
# tau, fold, binomial
t1m = list(TAUS.t1); t1m[5] = t1m[5] + GQ(Fr(1,100))
run('M19 tau1_5 perturbed in candidate ehat only', cand_kw={'taus': Taus(t1m, TAUS.t2)}, t1_kw={'taus_for_T1': TAUS})
t2m = list(TAUS.t2); t2m[3] = t2m[3] + GQ(0, Fr(1,100))
run('M20 tau2_3 perturbed in T1 operator only', t1_kw={'taus_for_T1': Taus(TAUS.t1, t2m)})
run('M21 fold sign flipped', t1_kw={'fold_sign': -1})
run('M22 fold binomial comb(mu-2,mu-m) -> comb(mu-1,mu-m)', t1_kw={'binom_fn': lambda a, b: comb(a+1, b)})
# E3 / diagonal initial
def shift(p, delta):
    return QSI(GQ(0), GQ(0)) if delta != 6 else E00(p, p).mul(GQ(Fr(1, 7)))
run('M23 diag INITIAL rho_pp(6) perturbed by rho_pp(0)/7 (E3 must catch via adjacent K=0)', cand_kw={'diag_initial_shift': shift})
def shift5(p, delta):
    return QSI(GQ(0), GQ(0)) if delta != 5 else E00(p, p).mul(GQ(0, Fr(1, 3)))
run('M24 diag INITIAL rho_pp(5) perturbed', cand_kw={'diag_initial_shift': shift5})
# E3 control: replacing the initial by the stable extrapolation sum_s (-1)^s ehat(s) must leave every adjacent cell passing
class CandExtrap(Cand):
    def get(self, p, q, k, delta):
        if p == q and k == 2*p:
            tot = ZERO
            for s in range(delta//2+1):
                tot = tot + self.Ehat(p, delta, s).mul(GQ((-1)**s))
            return tot
        return Cand.get(self, p, q, k, delta)
c = CandExtrap(dict(base))
tested, bad, per = t1_test(c, DL)
print(f'[E3 CONTROL] diag initial := stable extrapolation sum (-1)^s ehat: cells {tested}, bad {len(bad)} {classify(bad)} -> blindness: children unaffected', flush=True)
# Delta_pp is genuinely nonzero (so E3 is not vacuous)
c0 = Cand(dict(base)); ce = CandExtrap(dict(base))
nz = []
for p in (3, 4, 5, 6, 8):
    for delta in range(0, 7):
        if delta > 4*p-4: continue
        diff = c0.get(p, p, 2*p, delta) - ce.get(p, p, 2*p, delta)
        if not diff.iszero(): nz.append((p, delta))
print('Delta_pp(delta) != 0 at (p,delta):', nz[:20], '... total', len(nz), 'of', 5*7)
