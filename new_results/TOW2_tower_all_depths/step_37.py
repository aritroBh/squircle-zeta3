import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)
from towlib import *
import sympy as sp

d, m = sp.symbols('d m')
i = sp.I
aq = sp.sqrt(2)+1; bq = sp.sqrt(2)-1

# t(v) derivatives at v=i: t(i) = t*, t'(i)=0, t^{(r)}(i) = (-1)^{r-1} r!/(4 i^{r+1}) for r>=2
def tder(r):
    if r == 1: return sp.Integer(0)
    return (-1)**(r-1) * sp.factorial(r) / (4 * i**(r+1))

# S log-derivative ratios R_l = S^{(l)}/S at t*
N = {1: 4*m + (-2+6*i)*d - 4*i*d**2,
     2: (-40*i)*m + 16*m**2 + (sp.Rational(80,3)+20*i)*d + (-16+96*i)*d*m + (-80-48*i)*d**2 + (-32*i)*d**2*m + (sp.Rational(208,3)+16*i)*d**3 - 16*d**4}
def Pi(l):
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*d-2*jj-1)
    return out
R = {l: N[l]/Pi(l) for l in N}

# unnormalized log derivatives of S at t*
sl = {1: R[1], 2: sp.expand(R[2] - R[1]**2)}

# partial Bell polynomial B_{r,l}(x_2, x_3, ...) with blocks of size >= 2 (x_1 = 0)
from sympy.functions.combinatorial.numbers import bell as bell_sym
def partial_bell_no1(r, l, xargs):
    # xargs: dict block size -> value (sizes >= 2)
    tot = sp.Integer(0)
    # partitions of r into l parts each >= 2
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
        # multiplicity factor r! / (prod block! * prod mult_j!)
        from collections import Counter
        C = Counter(part)
        denom = sp.Integer(1)
        for sz in part: denom *= sp.factorial(sz)
        for sz, mult in C.items(): denom *= sp.factorial(mult)
        tot += sp.factorial(r)/denom * sp.prod(xargs[sz] for sz in part)
    return tot

xargs = {r: tder(r) for r in range(2, 8)}

def Lr(r):
    """(log A)^{(r)}(i)/r! as rational function of (d,m)."""
    tot = (2*d-1-0)*( (-1)**(r-1) * i**(-r) / r )   # (2d-1-D) log v part... note D handled separately below
    return tot

# Build L_r properly: pieces
def Lr_full(r, Dsym):
    # (2d-1-D) log v: (2d-1-D)(-1)^{r-1} i^{-r}/r
    p1 = (2*d-1-Dsym) * (-1)**(r-1) * i**(-r) / r
    # D log(4v): D (-1)^{r-1} i^{-r}/r
    p2 = Dsym * (-1)**(r-1) * i**(-r) / r
    # logS(t(v)): Bell
    p3 = sp.Integer(0)
    for l in range(1, r//2+1):
        p3 += sl[l] * partial_bell_no1(r, l, xargs) / sp.factorial(r)
    # (2-2d) log(v+i): (2-2d)(-1)^{r-1}(2i)^{-r}/r
    p4 = (2-2*d) * (-1)**(r-1) * (2*i)**(-r) / r
    # -log(v^2-a^2)-log(v^2-b^2)
    s5 = ( (i-aq)**(-r) + (i+aq)**(-r) + (i-bq)**(-r) + (i+bq)**(-r) )
    p5 = -(-1)**(r-1) * s5 / r
    return sp.expand(p1+p2+p3+p4+p5)

# D = deg S_m: from data, D = 2m? no... deg S_m(t) as polynomial in t. Let me check: S_poly(d,m) degree.
from genintegrand import S_poly
for (dd,mm) in [(5,5),(5,3),(6,4),(7,2)]:
    P = S_poly(dd,mm)
    print(f'deg S_{mm} (d={dd}):', len(P)-1)
