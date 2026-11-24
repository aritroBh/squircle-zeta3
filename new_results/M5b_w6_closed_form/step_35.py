import sys, os
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from flint import acb, arb, ctx
import m5b_precision as MP
from math import comb, factorial
# replicate run('w6',600) manually
ctx.prec = 600
eps_exp = 700; eps = arb(2)**(-eps_exp)
r = acb(2).sqrt(); a, b = 1+r, r-1
w = eps_exp*arb(2).log()
import json
d = json.load(open('w6_600_partial.json'))
pieces = {k: acb(arb(v[0], arb(v[1])), arb(v[2], arb(v[3]))) for k,v in d.items()}
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
print('Cb Mb Ca La finite:', Cb.is_finite(), Mb.is_finite(), Ca.is_finite(), La.is_finite())
T = [(i, j, [(z, acb(str(cp)) + r*acb(str(cq))) for z, cp, cq in ts]) for i, j, ts in MP.load_terms('w6')]
poles = {'0': acb(0), 'a': a, '-a': -a, 'b': b, '-b': -b}
tail_b, tail_a = arb(0), arb(0)
for i, j, ts in T:
    pb = sum((mag(c/(vb-poles[z])) for z, c in ts), arb(0))
    pa = sum((mag(c/(va-poles[z])) for z, c in ts), arb(0))
    tail_b += pb*moment(Cb, i)*Mb**j
    tail_a += pa*La**i*moment(Ca, j)
tail = (tail_a+tail_b).abs_upper()
print('tail:', tail.str(6), tail.is_finite())
value = sum(pieces.values(), acb(0))+acb(arb(0, tail))
print('value finite:', value.is_finite(), 'imag contains 0:', value.imag.contains(0))
print(value.real.rad().str(8))
