# step_64: redo PSLQ for R0 with OCT9c's settings; resolve inv0 zeta(3)
import os, sys, pickle
from fractions import Fraction as Fr
import mpmath as mp
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT9b_octagon_debug"))
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT9c_octagon_close"))
import oct9c_common as oc
import sympy as sp
D = oc.load_all()
engine = oc.get_engine(0)
mp.mp.dps = 250
av = oc.atom_values(engine, 0, D['lam'], 250)
B0 = D['brackets']['B0']
sym2key = engine.ATOMS
W3e, Wloe = oc.split_weight(B0, sym2key)
def ev(expr):
    tot = mp.mpf(0)
    for term in sp.Add.make_args(sp.expand(expr)):
        coeff, rest = term.as_coeff_Mul()
        v = mp.mpf(coeff.numerator)/mp.mpf(coeff.denominator) if coeff.is_rational else mp.mpf(str(coeff))
        for f in sp.Mul.make_args(rest):
            b_, e_ = f.as_base_exp()
            v *= av[b_]**int(e_)
        tot += v
    return tot
W3v = ev(W3e)
s0 = pickle.load(open(os.path.join(base,"research/08_theory/OCT9b_octagon_debug/sol_B0.pkl"),"rb"))
sols, negw7, logprods7 = oc.rhs_terms_2357()
solB0 = [sp.Rational(s) for s in s0['solB0']]
r_rhs = oc.eval_rhs_2357(solB0, sols, negw7, logprods7, 230)
R0 = W3v - r_rhs
print("R0 =", mp.nstr(R0, 50))
consts = [mp.zeta(3)] + [mp.pi**2/6*mp.log(p) for p in (2,3,5,7)]
rel = mp.pslq([R0]+consts, tol=mp.mpf('1e-210'), maxcoeff=10**20, maxsteps=10000)
print("pslq:", rel)
if rel:
    from fractions import Fraction as Fr
    c0 = Fr(rel[0]); parts = [Fr(c)/c0 for c in rel[1:]]
    print("R0 = -(", parts, ") . consts")
    print("zeta3 coefficient:", -parts[0], " vs OCT9c inv0:", Fr(30167467635454913,654492976038000))
