import sys, json, os
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from flint import acb, arb, ctx
from math import comb, factorial
ctx.prec = 600
eps_exp = 700
eps = arb(2)**(-eps_exp)
r = acb(2).sqrt(); a, b = 1+r, r-1
w = eps_exp*arb(2).log()
def mag(z): return abs(z.real)+abs(z.imag)
vb = b+acb(arb(0, eps)); va = a+acb(arb(0, eps))
import m5b_precision as MP
T = [(i, j, [(z, acb(str(cp)) + r*acb(str(cq))) for z, cp, cq in ts]) for i, j, ts in MP.load_terms('w6')]
poles = {'0': acb(0), 'a': a, '-a': -a, 'b': b, '-b': -b}
for i, j, ts in T:
    for z, c in ts:
        mm = mag(c/(vb-poles[z]))
        print(i, j, z, mm.str(6), mm.is_finite())
