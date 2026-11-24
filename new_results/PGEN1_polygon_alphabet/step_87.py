# step_59: recompute OCT9c edge-0 invisible zeta(3) from stored data
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
mp.mp.dps = 120
av = oc.atom_values(engine, 0, D['lam'], 120)
W3e, Wloe = D['Wp0']
# evaluate W3_e0 and Wlo_e0
sym2key = engine.ATOMS
def ev(expr):
    expr = sp.expand(expr)
    tot = mp.mpf(0)
    for term in sp.Add.make_args(expr):
        coeff, rest = term.as_coeff_Mul()
        v = mp.mpf(float(sp.Rational(coeff))) if coeff.is_rational else mp.mpf(str(coeff))
        for f in sp.Mul.make_args(rest):
            b_, e_ = f.as_base_exp()
            v *= av[b_]**int(e_)
        tot += v
    return tot
W3v = ev(W3e)
print("W3_e0 =", mp.nstr(W3v, 30))
print("Wlo_e0 =", mp.nstr(ev(Wloe), 30))
# Li3 part (OCT9c solB0), products, logprods
s0 = pickle.load(open(os.path.join(base,"research/08_theory/OCT9b_octagon_debug/sol_B0.pkl"),"rb"))
sols, negw7, logprods7 = oc.rhs_terms_2357()
solB0 = [sp.Rational(s) for s in s0['solB0']]
r_rhs = oc.eval_rhs_2357(solB0, sols, negw7, logprods7, 100)
print("OCT9c RHS0 (Li3+prod+log3) =", mp.nstr(r_rhs, 30))
R0 = W3v - r_rhs
print("R0 =", mp.nstr(R0, 30))
print("OCT9c quoted inv0 zeta3-part  =", mp.nstr(mp.mpf(30167467635454913)/mp.mpf(654492976038000)*mp.zeta(3),30))
