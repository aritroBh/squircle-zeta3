# step_78: edge-0 L3 via CR presentation vs OCT9c vs OCT9d
import pickle
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 120
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(q):
    q = mpfq(q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(abs(1-q))
z3 = mp.zeta(3)
CRp = pickle.load(open("CR_presentations.pkl","rb"))
pres = CRp[0]['presentation']
S_CR = mp.fsum([mpfq(c)*L3(z) for c,z in pres])
print("S_CR =", mp.nstr(S_CR, 40))
import os
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
Sc = mp.fsum([mpfq(c)*L3(z) for c,z in zip(solB0[:len(sols)], sols) if c])
B0a = [Fr(s) for s in sf['B0']]
sol0 = [Fr(s) for s in sf['sol0']]
Sd = mp.fsum([mpfq(c)*L3(z) for c,z in zip(sol0[:20], B0a) if c])
print("Sc   =", mp.nstr(Sc, 40))
print("Sd   =", mp.nstr(Sd, 40))
print("S_CR - Sd =", mp.nstr(S_CR-Sd, 15), " -> /z3:", mp.nstr((S_CR-Sd)/z3, 25))
print("S_CR - Sc =", mp.nstr(S_CR-Sc, 15), " -> /z3:", mp.nstr((S_CR-Sc)/z3, 25))
import sympy as sp
print("nsimplify (S_CR-Sd)/z3:", sp.nsimplify(str((S_CR-Sd)/z3)))
print("nsimplify (S_CR-Sc)/z3:", sp.nsimplify(str((S_CR-Sc)/z3)))
