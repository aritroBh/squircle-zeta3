# step_69: is (Sc - Sd)/zeta3 rational? verify symbols of both edge-0 presentations directly
import pickle, os, sys
from fractions import Fraction as Fr
import mpmath as mp
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
for dps in (60, 120):
    mp.mp.dps = dps
    def Li3(x): return mp.re(mp.polylog(3, x))
    def Li2(x): return mp.re(mp.polylog(2, x))
    def L3mpf(q):
        if q > 1: return L3mpf(1/q)
        return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(1-q)
    sols = [Fr(s) for s in s0['meta']['sols']]
    solB0 = [Fr(s) for s in s0['solB0']]
    Sc = mp.fsum([mpfq(c)*L3mpf(mpfq(z)) for c,z in zip(solB0[:len(sols)], sols) if c])
    B0a = [Fr(s) for s in sf['B0']]
    sol0 = [Fr(s) for s in sf['sol0']]
    Sd = mp.fsum([mpfq(c)*L3mpf(mpfq(z)) for c,z in zip(sol0[:20], B0a) if c])
    q = (Sc-Sd)/mp.zeta(3)
    print(f"dps={dps}: Sc-Sd = {mp.nstr(Sc-Sd, 30)}")
    print(f"   q = (Sc-Sd)/z3 = {mp.nstr(q, 40)}")
import sympy as sp
mp.mp.dps = 150
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3mpf(q):
    if q > 1: return L3mpf(1/q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(1-q)
sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
Sc = mp.fsum([mpfq(c)*L3mpf(mpfq(z)) for c,z in zip(solB0[:len(sols)], sols) if c])
B0a = [Fr(s) for s in sf['B0']]
sol0 = [Fr(s) for s in sf['sol0']]
Sd = mp.fsum([mpfq(c)*L3mpf(mpfq(z)) for c,z in zip(sol0[:20], B0a) if c])
q = (Sc-Sd)/mp.zeta(3)
print("nsimplify q:", sp.nsimplify(str(q)))
