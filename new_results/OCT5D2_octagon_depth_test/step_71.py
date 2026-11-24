import os, sys, pickle
os.chdir("<project>/research/08_theory/OCT5D2_octagon_depth_test")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import symengine as S
from symengine import Q, ZERO, ONE, add, scale, gsym
from fractions import Fraction as F
R2 = Q(0,1); B = Q(-1,1)
def w32(xx, yy): return (ZERO, ZERO, ONE/xx, ZERO, ONE/(xx*yy))
LI32 = [('(b,1)',(B,ONE)),('(1/sqrt2,1)',(R2/2,ONE)),('(-b,1)',(-B,ONE)),
        ('(b,-1)',(B,-ONE)),('(-b,-1)',(-B,-ONE)),('(1/sqrt2,-1)',(R2/2,-ONE)),
        ('(1-1/sqrt2,1)',(ONE-R2/2,ONE)),('(2-sqrt2,1)',(2-R2,ONE)),('(-b^2,1)',(-B*B,ONE))]
syms = []
words = set()
for name,(x,y) in LI32:
    s = gsym(w32(x,y))
    syms.append(s)
    words |= set(s)
words = sorted(words)
# matrix atoms x words, find kernel over Q (coeffs rational — atom symbols have rational coeffs? 
# check: are all symbol coefficients rational (q=0)? 
print("all Li32 symbol coeffs rational:", all(c.q==0 for s in syms for c in s.values()))
Mm = [[s.get(w, ZERO).p for s in syms] for w in words]  # rows words, cols atoms
import sympy as sp
Ms = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in Mm])
ns = Ms.nullspace()
print("symbol rank of 9 Li32 atoms:", Ms.rank(), "kernel dim:", len(ns))
for vec in ns:
    combo = [sp.nsimplify(v) for v in vec]
    print("kernel combo:", [(LI32[i][0], str(combo[i])) for i in range(9) if combo[i] != 0])
