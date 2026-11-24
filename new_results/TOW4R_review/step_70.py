
import os, sys, pickle, time
import sympy as sp
from fractions import Fraction as Fr
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
d, m = sp.symbols('d m'); i = sp.I
ef = pickle.load(open(work + "/tow4_copy/e_fits.pkl", 'rb'))
rec = ef[(7,3)]
mons, sol = rec['mons'], rec['sol']
poly = sum((re+im*i)*d**dm*m**mm for (dm,mm),(re,im) in zip(mons, sol))
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out
target = -2*i*(d-5)/Pi(3)
print("stale fit (7,3) (b=0, polynomial):", sp.expand(poly))
print("fit == -2i(d-5)/Pi3 ?", sp.expand(poly - target) == 0)
print("fit == 5i(d-5)/Pi3 (truth) ?", sp.expand(poly - 5*i*(d-5)/Pi(3)) == 0)
