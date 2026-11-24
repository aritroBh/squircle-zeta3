
import os, sys, pickle
import sympy as sp
from fractions import Fraction as Fr
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
d, m = sp.symbols('d m'); i = sp.I
ef = pickle.load(open(work + "/tow4_copy/e_fits.pkl", 'rb'))
rec = ef[(7,3)]
mons = eval(rec['mons']); sol = eval(rec['sol'])
print("b =", rec['b'], "deg =", rec['deg'], "npts =", rec['npts'])
poly = sum((re+im*i)*d**dm*m**mm for (dm,mm),(re,im) in zip(mons, sol))
b = int(rec['b'])
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out
cand = sp.expand(poly/Pi(b)) if b>0 else sp.expand(poly)
print("fit (7,3):", sp.cancel(cand) if b>0 else cand)
target = -2*i*(d-5)/Pi(3)
print("fit == -2i(d-5)/Pi3 (the misdisplayed value)?", sp.expand(cand - target) == 0)
# zigzag truth:
print("proved F(7,3) = 5i(d-5)/Pi3")
