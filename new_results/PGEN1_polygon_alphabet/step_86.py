# step_58: task 3 computation (mpf-safe recursion)
import pickle, os
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 100
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s123 = pickle.load(open(os.path.join(O9B,"sol_B123.pkl"),"rb"))
lists = pickle.load(open(os.path.join(O9B,"Sunit_lists.pkl"),"rb"))
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
def mpfq(q): q = Fr(q); return mp.mpf(q.numerator)/mp.mpf(q.denominator)
def Li3(x): return mp.re(mp.polylog(3, x))
def Li2(x): return mp.re(mp.polylog(2, x))
def L3mpf(q):
    """q: mpf real, q != 0,1. Zagier L3: direct for q<1, inversion L3(x)=L3(1/x) for q>1."""
    if q > 1: return L3mpf(1/q)
    return Li3(q) - mp.log(abs(q))*Li2(q) - mp.mpf(1)/3*mp.log(abs(q))**2*mp.log(1-q)
def L3rat(q): return L3mpf(mpfq(q))
z3 = mp.zeta(3)
Gargs = [Fr(a) for a in lists['Gargs']]
sols0 = [Fr(s) for s in s0['meta']['sols']]
G8 = [Fr(s) for s in sf['G8']]
B0 = [Fr(s) for s in sf['B0']]
invc = {0: Fr(30167467635454913,654492976038000), 1: Fr(113678924735323,3700404000000),
        2: Fr(1978761997948319,41624512450560), 3: Fr(15326623217074655,142615945728000)}
invd = {0: Fr(-7340,5043), 1: Fr(34497,15625), 2: Fr(2085,2197), 3: Fr(-57300,24389)}
results = {}
for e in range(4):
    if e == 0:
        Xc = mp.fsum([mpfq(c)*L3rat(z) for c,z in zip([Fr(s) for s in s0['solB0'][:len(sols0)]], sols0) if c]) + mpfq(invc[0])*z3
        Xd = mp.fsum([mpfq(c)*L3rat(z) for c,z in zip([Fr(s) for s in sf['sol0'][:20]], B0) if c]) + mpfq(invd[0])*z3
    else:
        Xc = mp.fsum([mpfq(c)*L3rat(z) for c,z in zip([Fr(s) for s in s123[e]['sol'][:len(Gargs)]], Gargs) if c]) + mpfq(invc[e])*z3
        Xd = mp.fsum([mpfq(c)*L3rat(z) for c,z in zip([Fr(s) for s in sf['sols123'][e][:8]], G8) if c]) + mpfq(invd[e])*z3
    results[e] = (Xc, Xd)
    print(f"edge {e}: L3 OCT9c = {mp.nstr(Xc, 40)}")
    print(f"        L3 OCT9d = {mp.nstr(Xd, 40)}")
    print(f"        |diff| = {mp.nstr(abs(Xc-Xd), 5)}")
pickle.dump(results, open("task3_L3.pkl","wb"))
