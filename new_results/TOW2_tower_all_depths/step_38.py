import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
import sympy as sp
from collections import Counter

d, m = sp.symbols('d m')
i = sp.I
aq = sp.sqrt(2)+1; bq = sp.sqrt(2)-1

def tder(r):
    if r == 1: return sp.Integer(0)
    return (-1)**(r-1) * sp.factorial(r) / (4 * i**(r+1))

N = {1: 4*m + (-2+6*i)*d - 4*i*d**2,
     2: (-40*i)*m + 16*m**2 + (sp.Rational(80,3)+20*i)*d + (-16+96*i)*d*m + (-80-48*i)*d**2 + (-32*i)*d**2*m + (sp.Rational(208,3)+16*i)*d**3 - 16*d**4}
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out
R = {l: N[l]/Pi(l) for l in N}
sl = {1: R[1], 2: sp.expand(R[2] - R[1]**2)}

def partial_bell_no1(r, l, xargs):
    tot = sp.Integer(0)
    def parts(rem, minpart, cur):
        if len(cur) == l:
            if rem == 0: yield tuple(sorted(cur))
            return
        for p in range(minpart, rem+1):
            yield from parts(rem-p, p, cur+[p])
    seen = set()
    for part in parts(r, 2, []):
        if part in seen: continue
        seen.add(part)
        C = Counter(part)
        denom = sp.Integer(1)
        for sz in part: denom *= sp.factorial(sz)
        for sz, mult in C.items(): denom *= sp.factorial(mult)
        tot += sp.factorial(r)/denom * sp.prod(xargs[sz] for sz in part)
    return tot

xargs = {r: tder(r) for r in range(2, 8)}
def Lr_full(r):
    Dsym = 2*d
    p12 = (2*d-1) * (-1)**(r-1) * i**(-r) / r    # (2d-1-D)log v + D log(4v) with D=2d
    p3 = sp.Integer(0)
    for l in range(1, r//2+1):
        p3 += sl[l] * partial_bell_no1(r, l, xargs) / sp.factorial(r)
    p4 = (2-2*d) * (-1)**(r-1) * (2*i)**(-r) / r
    s5 = ( (i-aq)**(-r) + (i+aq)**(-r) + (i-bq)**(-r) + (i+bq)**(-r) )
    p5 = -(-1)**(r-1) * s5 / r
    return sp.expand(p12+p3+p4+p5)

L = {r: Lr_full(r) for r in range(1,5)}
print('L1 =', sp.simplify(L[1]), ' (should be -i(d-2) =', -i*(d-2), ')')
for r in range(2,5):
    print(f'L{r} =', sp.simplify(L[r]))

# complete Bell polynomials for rho(delta)/rho(0) = [y^delta] exp(sum L_r y^r)
def complete_bell(delta, Ls):
    tot = sp.Integer(0)
    def parts(rem, minpart, cur):
        if rem == 0:
            yield tuple(sorted(cur)); return
        for p in range(minpart, rem+1):
            yield from parts(rem-p, p, cur+[p])
    for part in parts(delta, 1, []):
        C = Counter(part)
        coef = sp.Integer(1)
        for sz, mult in C.items():
            coef *= (Ls[sz]**mult) / sp.factorial(mult)
        tot += coef
    return tot

a = {delta: sp.expand(complete_bell(delta, L)) for delta in range(1,5)}
for delta in range(1,5):
    print(f'a_{delta} = rho({delta})/rho(0) =', sp.simplify(a[delta]))
import pickle
pickle.dump({'L': {r: str(L[r]) for r in L}, 'a': {de: str(a[de]) for de in a}}, open('rho_closed_forms.pkl','wb'))
