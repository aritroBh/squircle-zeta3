import sys, os
sys.dont_write_bytecode=True
sys.path.insert(0, os.getcwd())
import sympy as sp
import prove_E3b78 as PE
T4=PE.T4
# negative control 1: perturb tau2(8) (enters the direct term of E3(8)) -> must NOT vanish
t1,t2=T4._TAU8[8]
T4._TAU8[8]=(t1, t2+sp.Rational(1,128))
print('E3(8) with perturbed tau2(8): zero?', PE.check_zero(PE.feed_terms(8)))
T4._TAU8[8]=(t1,t2)
# negative control 2: perturb the discrepancy Delta(3) by a constant -> the feed at delta=7 must not vanish
PE._D.clear(); PE._Dc.clear(); PE._g.clear()
orig=PE.Delta
def bad_delta(delta, conj=False):
    v=orig(delta,conj)
    return v + (sp.Rational(1,7) if delta==3 else 0)
PE.Delta=bad_delta
PE._Dc.clear(); PE._g.clear()
print('E3(7) with Delta(3) perturbed: zero?', PE.check_zero(PE.feed_terms(7)))
