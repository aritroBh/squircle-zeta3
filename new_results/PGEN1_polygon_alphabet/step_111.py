# step_83: direct check: is sum_c s(z_c) - sum_d s(z_d) in D? plus series L3 cross-check
import os, sys, pickle
from fractions import Fraction as Fr
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT8_symbol_test"))
from oct8 import Quotient, S_Li3, vec_of_rational as vf
import sympy as sp
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
sf = pickle.load(open(os.path.join(O9D,"short_form.pkl"),"rb"))
Q = Quotient((2,3,5,7))
def coord(S):
    v = sp.zeros(Q.dim3,1)
    for t,c in S.items():
        v[Q.bidx[t],0] += sp.Rational(c.numerator, c.denominator)
    return v
sols = [Fr(s) for s in s0['meta']['sols']]
solB0 = [Fr(s) for s in s0['solB0']]
B0a = [Fr(s) for s in sf['B0']]
sol0 = [Fr(s) for s in sf['sol0']]
Sc_col = sp.zeros(Q.dim3,1)
for c,z in zip(solB0[:len(sols)], sols):
    if c: Sc_col += sp.Rational(c.numerator,c.denominator)*coord(S_Li3(z, vf))
Sd_col = sp.zeros(Q.dim3,1)
for c,z in zip(sol0[:20], B0a):
    if c: Sd_col += sp.Rational(c.numerator,c.denominator)*coord(S_Li3(z, vf))
diff = Sc_col - Sd_col
MD = sp.Matrix.hstack(*[d.T for d in Q.D])
print("Li3-parts difference in D:", sp.Matrix.hstack(MD, diff).rank() == MD.rank(), " rank D =", MD.rank(), " rank with diff:", sp.Matrix.hstack(MD, diff).rank())
# series cross-check of L3 at a few arguments
import mpmath as mp
mp.mp.dps = 100
def Li3_series(z, terms=200000):
    z = mp.mpf(z)
    if abs(z) >= 1: return None
    return mp.nsum(lambda n: z**int(n)/mp.mpf(int(n))**3, [1, mp.inf])
def L3_series(z):
    z = mp.mpf(z)
    li3 = Li3_series(z)
    li2 = mp.nsum(lambda n: z**int(n)/mp.mpf(int(n))**2, [1, mp.inf])
    return li3 - mp.log(abs(z))*li2 - mp.mpf(1)/3*mp.log(abs(z))**2*mp.log(1-z)
def L3_pl(z):
    z = mp.mpf(z)
    return mp.re(mp.polylog(3,z)) - mp.log(abs(z))*mp.re(mp.polylog(2,z)) - mp.mpf(1)/3*mp.log(abs(z))**2*mp.log(abs(1-z))
for z in ['0.1','-0.4','0.75','-0.9','0.5','-0.05']:
    a, b = L3_series(z), L3_pl(z)
    print(f"L3({z}): series={mp.nstr(a,25)} polylog={mp.nstr(b,25)} diff={mp.nstr(abs(a-b),4)}")
