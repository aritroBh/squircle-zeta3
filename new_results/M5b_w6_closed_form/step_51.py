import sys, os, json
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from flint import acb, arb, ctx
import m5b_precision as MP
ctx.prec = 4800
bits=4800
eps_exp = bits+100
r = acb(2).sqrt(); a, b = 1+r, r-1
w = eps_exp*arb(2).log()
T = [(i, j, [(z, acb(str(cp)) + r*acb(str(cq))) for z, cp, cq in ts]) for i, j, ts in MP.load_terms('w5pp')]
poles = {'0': acb(0), 'a': a, '-a': -a, 'b': b, '-b': -b}
def integrand(x, analytic=False, side='left'):
    e = x.exp(); v = b+e
    l = x+((v+a)/(4*v)).log(analytic=analytic)
    m = ((a-v)*(v+b)/(4*v)).log(analytic=analytic)
    tot = acb(0)
    for i, j, ts in T:
        omv = acb(0)
        for z, c in ts:
            omv += c/(v-poles[z])
        tot += (l**i)*(m**j)*omv
    return e*tot
lo, hi = -acb(w), (1-b).log()
print('interval', lo.str(8), hi.str(8))
v1 = acb.integral(lambda x, analytic: integrand(x, analytic), lo, hi)
print('one-shot:', v1.real.mid().str(30), 'rad', v1.real.rad().str(6))
mid = acb(-40)
v2 = acb.integral(lambda x, analytic: integrand(x, analytic), lo, mid) + \
     acb.integral(lambda x, analytic: integrand(x, analytic), mid, hi)
print('split@-40:', v2.real.mid().str(30), 'rad', v2.real.rad().str(6))
