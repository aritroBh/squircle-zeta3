# step_86: pairwise L3 differences with proper pslq at 200 digits
import pickle, os
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 200
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
CRp = pickle.load(open("CR_presentations.pkl","rb"))
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(q):
    q = mpfq(q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(abs(1-q))
z3 = mp.zeta(3)
sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
Sc = mp.fsum([mpfq(c)*L3(z) for c,z in zip(solB0[:len(sols)], sols) if c])
B0a = [Fr(s) for s in sf['B0']]
sol0 = [Fr(s) for s in sf['sol0']]
Sd = mp.fsum([mpfq(c)*L3(z) for c,z in zip(sol0[:20], B0a) if c])
S_CR = mp.fsum([mpfq(c)*L3(z) for c,z in CRp[0]['presentation']])
for nm, val in [("Sc-Sd", Sc-Sd), ("S_CR-Sd", S_CR-Sd), ("S_CR-Sc", S_CR-Sc)]:
    q = val/z3
    rel = mp.pslq([q, mp.mpf(1)], tol=mp.mpf('1e-150'), maxcoeff=10**35, maxsteps=10000)
    print(f"{nm}/z3 = {mp.nstr(q, 40)}")
    print(f"   pslq: {rel}")
    if rel:
        print(f"   rational: {Fr(rel[1])/Fr(rel[0])} = {Fr(rel[1])/Fr(rel[0])*(-1)}")
