# step_66: edge 0: number identity check Li3-side for OCT9c vs OCT9d
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
mp.mp.dps = 100
av = oc.atom_values(engine, 0, D['lam'], 100)
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
print("W3_e0 =", mp.nstr(W3v, 25))
# OCT9c: r_rhs = Lam+Pi+log3 ; R0 = inv0 (zeta3 + zeta2logp)
s0 = pickle.load(open(os.path.join(base,"research/08_theory/OCT9b_octagon_debug/sol_B0.pkl"),"rb"))
sols, negw7, logprods7 = oc.rhs_terms_2357()
solB0 = [sp.Rational(s) for s in s0['solB0']]
r_rhs_c = oc.eval_rhs_2357(solB0, sols, negw7, logprods7, 100)
invc0 = Fr(30167467635454913,654492976038000)
z3 = mp.zeta(3)
Xc_total = r_rhs_c + mp.mpf(invc0.numerator)/mp.mpf(invc0.denominator)*z3
print("OCT9c: rhs + inv0_z3 =", mp.nstr(Xc_total, 25))
print("       (should be W3_e0 - inv0_zeta2logp parts)")
inv0_z2lp = (Fr(18544071891072555,654492976038000)*mp.pi**2/6*mp.log(2)
            + Fr(40489125837856272,654492976038000)*mp.pi**2/6*mp.log(3)
            + Fr(15850884944620440,654492976038000)*mp.pi**2/6*mp.log(5)
            - Fr(2646857348685000,654492976038000)*mp.pi**2/6*mp.log(7))
print("W3 - zeta2logp part =", mp.nstr(W3v - inv0_z2lp, 25))
print("diff =", mp.nstr(Xc_total - (W3v - inv0_z2lp), 8))
