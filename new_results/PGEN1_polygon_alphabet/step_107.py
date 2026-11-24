# step_79: edge-1 and edge-2 L3 sums: CR presentation vs OCT9d presentation
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
import os
base = "<project>"
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
G8 = [Fr(s) for s in sf['G8']]
import sympy as sp
for e in (1,2,3):
    S_CR = mp.fsum([mpfq(c)*L3(z) for c,z in CRp[e]['presentation']])
    S_d = mp.fsum([mpfq(c)*L3(z) for c,z in zip([Fr(s) for s in sf['sols123'][e][:8]], G8) if c])
    q = (S_CR - S_d)/z3
    print(f"edge {e}: S_CR-Sd = {mp.nstr(S_CR-S_d,12)}  q = {mp.nstr(q,25)}  nsimplify: {sp.nsimplify(str(q))}")
