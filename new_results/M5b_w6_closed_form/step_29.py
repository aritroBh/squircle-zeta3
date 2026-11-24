import sys, json
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
vb = b+acb(arb(0, eps))
Cb = mag(((vb+a)/(4*vb)).log())
Mb = mag(((a-vb)*(vb+b)/(4*vb)).log())
print('Cb', Cb.str(5), 'Mb', Mb.str(5))
def moment(C, k):
    out, power = arb(0), arb(1)
    for j in range(k+1):
        out += arb(comb(k, j)*factorial(k-j))*power
        power *= w+C
    return eps*out
mo = moment(Cb, 3)
print('moment', mo.str(5), mo.is_finite())
import m5b_precision as MP
T = [(i, j, [(z, acb(str(cp)) + r*acb(str(cq))) for z, cp, cq in ts]) for i, j, ts in MP.load_terms('w6')]
poles = {'0': acb(0), 'a': a, '-a': -a, 'b': b, '-b': -b}
tail_b = arb(0)
for i, j, ts in T:
    pb = sum((mag(c/(vb-poles[z])) for z, c in ts), arb(0))
    print(i,j,'pb', pb.str(5), pb.is_finite())
    term = pb*moment(Cb, i)*Mb**j
    print('  term', term.str(5), term.is_finite())
    tail_b += term
print('tail_b', tail_b.str(5))
