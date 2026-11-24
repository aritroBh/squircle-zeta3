import sys, os, json
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from flint import acb, arb, ctx
from math import comb, factorial
import m5b_precision as MP
ctx.prec = 600
eps_exp = 700; eps = arb(2)**(-eps_exp)
r = acb(2).sqrt(); a, b = 1+r, r-1
w = eps_exp*arb(2).log()
def mag(z): return abs(z.real)+abs(z.imag)
def moment(C,k):
    out, power = arb(0), arb(1)
    for j in range(k+1):
        out += arb(comb(k,j)*factorial(k-j))*power
        power *= w+C
    return eps*out
vb, va = b+acb(arb(0, eps)), a+acb(arb(0, eps))
Cb = mag(((vb+a)/(4*vb)).log()); Mb = mag(((a-vb)*(vb+b)/(4*vb)).log())
Ca = mag(((va+b)/(4*va)).log()); La = mag(((va-b)*(va+a)/(4*va)).log())
T = [(i, j, [(z, acb(str(cp)) + r*acb(str(cq))) for z, cp, cq in ts]) for i, j, ts in MP.load_terms('w6')]
poles = {'0': acb(0), 'a': a, '-a': -a, 'b': b, '-b': -b}
for i, j, ts in T:
    pb = sum((mag(c/(vb-poles[z])) for z, c in ts), arb(0))
    pa = sum((mag(c/(va-poles[z])) for z, c in ts), arb(0))
    tb = pb*moment(Cb, i)*Mb**j
    ta = pa*La**i*moment(Ca, j)
    print(i,j,'tb', tb.str(6), tb.is_finite(), 'ta', ta.str(6), ta.is_finite())
