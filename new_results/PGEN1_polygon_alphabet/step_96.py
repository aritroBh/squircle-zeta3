# step_68: task 3 edge 0 careful recompute at 150 dps
import pickle, os
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 150
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3mpf(q):
    if q > 1: return L3mpf(1/q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(1-q)
z3 = mp.zeta(3)
sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
nL = len(sols)
Sc = mp.fsum([mpfq(c)*L3mpf(mpfq(z)) for c,z in zip(solB0[:nL], sols) if c])
print("sum_c L3 =", mp.nstr(Sc, 40))
B0a = [Fr(s) for s in sf['B0']]
sol0 = [Fr(s) for s in sf['sol0']]
Sd = mp.fsum([mpfq(c)*L3mpf(mpfq(z)) for c,z in zip(sol0[:20], B0a) if c])
print("sum_d L3 =", mp.nstr(Sd, 40))
invc0 = Fr(30167467635454913,654492976038000)
invd0 = Fr(-7340,5043)
lhs = Sc + mpfq(invc0)*z3
rhs = Sd + mpfq(invd0)*z3
print("Xc =", mp.nstr(lhs, 40))
print("Xd =", mp.nstr(rhs, 40))
print("Xc - Xd =", mp.nstr(lhs-rhs, 10))
print("(Xc-Xd)/z3 =", mp.nstr((lhs-rhs)/z3, 20))
# how many nonzero coeffs in c presentation
print("nonzero c:", sum(1 for c in solB0[:nL] if c), "of", nL)
