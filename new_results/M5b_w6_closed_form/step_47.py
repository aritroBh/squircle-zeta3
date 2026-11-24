import sys, os, json
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from flint import acb, arb, ctx
from math import comb, factorial
import m5b_precision as MP
ctx.prec = 4600
eps_exp = 4700; eps = arb(2)**(-eps_exp)
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
print('Cb', Cb.is_finite(), 'Mb', Mb.is_finite(), 'Ca', Ca.is_finite(), 'La', La.is_finite())
print('Cb', Cb.str(6), 'Ca', Ca.str(6))
T = [(i, j, [(z, acb(str(cp)) + r*acb(str(cq))) for z, cp, cq in ts]) for i, j, ts in MP.load_terms('w5pp')]
poles = {'0': acb(0), 'a': a, '-a': -a, 'b': b, '-b': -b}
def pow_int(x,k):
    out=arb(1)
    for _ in range(k): out*=x
    return out
for i, j, ts in T:
    pb = sum((mag(c/(vb-poles[z])) for z, c in ts), arb(0))
    pa = sum((mag(c/(va-poles[z])) for z, c in ts), arb(0))
    tb = pb*moment(Cb, i)*pow_int(Mb, j)
    ta = pa*pow_int(La, i)*moment(Ca, j)
    print(i,j,'tb', tb.str(6), tb.is_finite(), 'ta', ta.str(6), ta.is_finite())
