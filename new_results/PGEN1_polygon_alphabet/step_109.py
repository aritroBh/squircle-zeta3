# step_81: rationality tests at 200 digits
import pickle, os
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 200
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3(q):
    q = mpfq(q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(abs(1-q))
sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
Sc = mp.fsum([mpfq(c)*L3(z) for c,z in zip(solB0[:len(sols)], sols) if c])
B0a = [Fr(s) for s in sf['B0']]
sol0 = [Fr(s) for s in sf['sol0']]
Sd = mp.fsum([mpfq(c)*L3(z) for c,z in zip(sol0[:20], B0a) if c])
z3 = mp.zeta(3)
q = (Sc-Sd)/z3
print("Sc-Sd =", mp.nstr(Sc-Sd, 50))
print("q =", mp.nstr(q, 60))
# proper rationality tests
rel = mp.pslq([q, mp.mpf(1)], tol=mp.mpf('1e-150'), maxcoeff=10**30, maxsteps=10000)
print("pslq [q,1]:", rel)
idq = mp.identify(q, ['zeta(3)'])
print("identify:", idq)
# Xc-Xd with invs
invc0 = Fr(30167467635454913,654492976038000); invd0 = Fr(-7340,5043)
Xc = Sc + mpfq(invc0)*z3; Xd = Sd + mpfq(invd0)*z3
print("Xc-Xd =", mp.nstr(Xc-Xd, 30))
r = (Xc-Xd)/z3
print("(Xc-Xd)/z3 =", mp.nstr(r, 60))
rel2 = mp.pslq([r, mp.mpf(1)], tol=mp.mpf('1e-150'), maxcoeff=10**30, maxsteps=10000)
print("pslq [r,1]:", rel2)
