# s15: TRANSLATION check for E3 (blindness).  prove_E3b56.feed(delta) is linear in Delta_pp.  Perturb Delta -> Delta + x (arbitrary
# Gaussian rationals x_delta') in the script AND replace the diagonal initial value by (rho_pp + e00_pp x) in T1 at K=0:
# the script's feed (units of e00_pp) must equal [T1(actual+x) - T1(stable extrapolation)] / e00_pp  for every x.
import sys, time, random
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import sympy as sp
from candlib import *
import prove_E3b56 as PE
random.seed(7)
base = base_exprs()
class CandExtrap(Cand):
    def get(self, p, q, k, delta):
        if p == q and k == 2*p:
            tot = ZERO
            for s in range(delta//2+1):
                tot = tot + self.Ehat(p, delta, s).mul(GQ((-1)**s))
            return tot
        return Cand.get(self, p, q, k, delta)
B = CandExtrap(dict(base))
pS = sp.Symbol('p')
orig_Delta = PE.Delta
def to_gq(z):
    z = sp.nsimplify(sp.expand(z)); re_, im_ = z.as_real_imag()
    return GQ(Fr(int(re_.p), int(re_.q)), Fr(int(im_.p), int(im_.q)))
ntest = 0; nbad = 0; nonzero = 0
for trial in range(3):
    X = {dl: GQ(Fr(random.randint(-9, 9), random.randint(1, 5)), Fr(random.randint(-9, 9), random.randint(1, 5))) for dl in range(0, 7)}
    Xs = {dl: sp.Rational(X[dl].a.numerator, X[dl].a.denominator) + sp.I*sp.Rational(X[dl].b.numerator, X[dl].b.denominator) for dl in X}
    def Delta_patched(delta, conj=False, _X=Xs):
        b = orig_Delta(delta, conj)
        return sp.cancel(b + (sp.conjugate(_X[delta]) if conj else _X[delta]))
    PE.Delta = Delta_patched
    A = Cand(dict(base), diag_initial_shift=lambda p, dl, _X=X: E00(p, p).mul(_X[dl]))
    for p in (3, 4, 5, 8):
        for delta in (5, 6):
            q = p-1
            lhs = T1rhs(A.get, A.taus, p, q, 2*p, delta) - T1rhs(B.get, B.taus, p, q, 2*p, delta)
            assert lhs.x.iszero()
            mine = lhs.y / E00(p, p).y
            script = to_gq(PE.feed(delta).subs(pS, p))
            ntest += 1
            if not mine.iszero(): nonzero += 1
            if mine != script:
                nbad += 1; print('MISMATCH', trial, p, delta, mine, script)
PE.Delta = orig_Delta
print(f'E3 translation: script feed(delta) vs T1-derived feed of (Delta+x): {ntest} comparisons, {nbad} mismatches; nonzero feeds (sensitivity) {nonzero}')
# and with x = 0 both vanish (the blindness itself), pointwise:
z = 0
for p in (3, 4, 5, 6, 9):
    for delta in (5, 6):
        z += (PE.feed(delta).subs(pS, p) != 0)
print('script feed(delta) at x=0, p in {3,4,5,6,9}, delta in {5,6}: nonzero count', z)
