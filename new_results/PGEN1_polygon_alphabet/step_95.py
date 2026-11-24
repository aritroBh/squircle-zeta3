# step_67: edge 0: OCT9d side number check
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
sf = pickle.load(open(os.path.join(base,"research/08_theory/OCT9d_octagon_short/short_form.pkl"),"rb"))
B0args = [sp.Rational(s) for s in sf['B0']]
PB0 = [(sp.Rational(w), int(p)) for w,p in sf['PB0']]
sol0 = [sp.Rational(s) for s in sf['sol0']]
def Li3v(a): a = mp.mpf(a.numerator)/mp.mpf(a.denominator); return mp.re(mp.polylog(3, a))
def Li2v(a): a = mp.mpf(a.numerator)/mp.mpf(a.denominator); return mp.re(mp.polylog(2, a))
nL = len(B0args)
Lam_d = mp.fsum([mp.mpf(c.numerator)/mp.mpf(c.denominator)*Li3v(z) for c,z in zip(sol0[:nL], B0args)])
Pi_d = mp.fsum([mp.mpf(c.numerator)/mp.mpf(c.denominator)*(-Li2v(w)*mp.log(p)) for c,(w,p) in zip(sol0[nL:nL+len(PB0)], PB0)])
print("Lam_d =", mp.nstr(Lam_d, 25))
print("Pi_d  =", mp.nstr(Pi_d, 25))
R_d = W3v - Lam_d - Pi_d
print("R_d = W3 - Lam - Pi =", mp.nstr(R_d, 25), " OCT9d says -19.4499")
# OCT9d relation: 25845375 R + 37607500 z3 + 495501375 z2l2 + 1398618 z2l3 - 48773505 z2l5 + 5972500 z2l7 = 0
rel = [Fr(25845375), Fr(37607500), Fr(495501375), Fr(1398618), Fr(-48773505), Fr(5972500)]
consts = [mp.zeta(3)] + [mp.pi**2/6*mp.log(p) for p in (2,3,5,7)]
val = mp.mpf(rel[0].numerator)/mp.mpf(rel[0].denominator)*R_d
for c, x in zip(rel[1:], consts):
    val += mp.mpf(c.numerator)/mp.mpf(c.denominator)*x
print("OCT9d relation value (should be 0):", mp.nstr(val, 8))
