import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
from oct5d3_basis import li32_symbol, B, H, R2, ONE, ZERO, Q, add, scale, obstruction
from fractions import Fraction as F

LI32 = [('(b,1)',(B,ONE)),('(1/sqrt2,1)',(H,ONE)),('(-b,1)',(-B,ONE)),
        ('(b,-1)',(B,-ONE)),('(-b,-1)',(-B,-ONE)),('(1/sqrt2,-1)',(H,-ONE)),
        ('(1-1/sqrt2,1)',(ONE-H,ONE)),('(2-sqrt2,1)',(2-R2,ONE)),('(-b^2,1)',(-B*B,ONE))]
syms = []
words = set()
for name,(x,y) in LI32:
    s = li32_symbol(x,y)
    syms.append(s); words |= set(s)
words = sorted(words)
print("all symbol coeffs rational:", all(c.q==0 for s in syms for c in s.values()))
import sympy as sp
Mm = [[sp.Rational(s.get(w, ZERO).p.numerator, s.get(w, ZERO).p.denominator) for s in syms] for w in words]
Ms = sp.Matrix(Mm)
ns = Ms.nullspace()
print("9 Li32 atoms: symbol rank", Ms.rank(), "kernel dim", len(ns))
for vec in ns:
    combo = [(LI32[i][0], sp.nsimplify(vec[i])) for i in range(9) if vec[i] != 0]
    print("  kernel:", combo)
# obstruction pair of each atom
for name,(x,y) in LI32:
    print("  pair", name, obstruction(li32_symbol(x,y)))
