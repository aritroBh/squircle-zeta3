import sys, os
work = '<project>/research/08_theory/TOW2_tower_all_depths'
sys.path.insert(0, work); os.chdir(work)

lib = '''"""tow2lib.py -- TOW2 library: closed forms for the full Newton e-family at all depths.
Contents:
  zigzag A_n; rho(delta) closed forms (Bell polynomials in log-derivatives of A at i,
    using the PROVED S-derivative lemmas S1,S3,S4 = R_0,R_1,R_2 and fitted R_3,R_4);
  zigzag top coefficients e(2s,s), e(2s+1,s);
  descending-recursion constructor e_desc(delta,s) (sympy rational functions of (d,m));
  Newton-feed operator Lambda(delta,s) in closed form.
Everything is a rational function of (d, m=first index) over Q(i, sqrt2-consts).
"""
import sympy as sp
from collections import Counter

d, m = sp.symbols('d m')
i = sp.I
aq = sp.sqrt(2)+1; bq = sp.sqrt(2)-1

# ---------------- zigzag numbers A_n: sech x + tanh = sum A_n x^n/n! (unsigned) ----------------
def zigzag(nmax):
    x = sp.symbols('x')
    ser = sp.series(sp.sech(x)+sp.tanh(x), x, 0, nmax+1).removeO().expand()
    return {n: abs(sp.Rational(ser.coeff(x, n))*sp.factorial(n)) for n in range(nmax+1)}
A = zigzag(20)

# ---------------- S-derivative ratios R_l = S^{(l)}(t*)/S(t*) ----------------
# R_1, R_2 PROVED for all (d,m) (prove_S12.py); R_3, R_4 fitted+verified (prove_S34.py pending)
N = {1: 4*m + (-2+6*i)*d - 4*i*d**2,
     2: (-40*i)*m + 16*m**2 + (sp.Rational(80,3)+20*i)*d + (-16+96*i)*d*m + (-80-48*i)*d**2 + (-32*i)*d**2*m + (sp.Rational(208,3)+16*i)*d**3 - 16*d**4,
     3: (-688)*m + (-672*i)*m**2 + 64*m**3 + (344-280*i)*d + (2768+672*i)*d*m + (-96+864*i)*d*m**2 + (-1384+1144*i)*d**2 + (-3264-864*i)*d**2*m + (-192*i)*d**2*m**2 + (1648-1984*i)*d**3 + (1408+192*i)*d**3*m + (-704+1600*i)*d**4 + (-192)*d**4*m + (96-544*i)*d**5 + (64*i)*d**6,
     4: 13824*i*m - 30784*m**2 - 6912*i*m**3 + 256*m**4 + (-sp.Rational(26576,5)-6912*i)*d + (30784-88320*i)*d*m + (72448+10368*i)*d*m**2 + (-512+6144*i)*d*m**3 + (sp.Rational(60592,3)+44160*i)*d**2 + (-72448+169856*i)*d**2*m + (-54528-9216*i)*d**2*m**2 + (-1024*i)*d**2*m**3 + (-48000-86656*i)*d**3 + (54784-148736*i)*d**3*m + (15872+1536*i)*d**3*m**2 + (sp.Rational(202624,3)+75904*i)*d**4 + (-15872+64768*i)*d**4*m + (-1536)*d**4*m**2 + (-sp.Rational(251904,5)-32640*i)*d**5 + (1536-13312*i)*d**5*m + (sp.Rational(57856,3)+6656*i)*d**6 + (1024*i)*d**6*m + (-3584-512*i)*d**7 + 256*d**8}
def Pi(l, dd=None):
    dd = dd or d
    out = sp.Integer(1)
    for jj in range(1, l+1): out *= (2*dd-2*jj-1)
    return out
R = {l: N[l]/Pi(l) for l in N}

# unnormalized log-derivatives of S at t*
sl = {1: R[1],
      2: sp.expand(R[2] - R[1]**2),
      3: sp.expand(R[3] - 3*R[1]*R[2] + 2*R[1]**3),
      4: sp.expand(R[4] - 4*R[1]*R[3] - 3*R[2]**2 + 12*R[1]**2*R[2] - 6*R[1]**4)}

# ---------------- log-derivatives L_r of A(v) = Phi(v)(v-i)^{2d-2} at v=i ----------------
def tder(r):
    if r == 1: return sp.Integer(0)
    return (-1)**(r-1) * sp.factorial(r) / (4 * i**(r+1))

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

_xargs = {r: tder(r) for r in range(2, 12)}
def Lr(r):
    p12 = (2*d-1) * (-1)**(r-1) * i**(-r) / r
    p3 = sp.Integer(0)
    for l in range(1, r//2+1):
        p3 += sl[l] * partial_bell_no1(r, l, _xargs) / sp.factorial(r)
    p4 = (2-2*d) * (-1)**(r-1) * (2*i)**(-r) / r
    s5 = ( (i-aq)**(-r) + (i+aq)**(-r) + (i-bq)**(-r) + (i+bq)**(-r) )
    p5 = -(-1)**(r-1) * s5 / r
    return sp.expand(p12+p3+p4+p5)

_L = {}
def L(r):
    if r not in _L: _L[r] = Lr(r)
    return _L[r]

def complete_bell(delta):
    Ls = {r: L(r) for r in range(1, delta+1)}
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

_a = {}
def arho(delta):
    """rho(delta)/rho(0) = a_delta (rational function of (d,m))."""
    if delta == 0: return sp.Integer(1)
    if delta not in _a: _a[delta] = sp.expand(complete_bell(delta))
    return _a[delta]

# ---------------- e00 ratios, tau constants ----------------
def rA0(dd=None, mm=None):
    """e00_{p+1,q}/e00_{p,q} = i(2d-1)/(p+1) (non-diagonal both)."""
    dd = dd or d; mm = mm or m
    return i*(2*dd-1)/(mm+1)
def rB0(dd=None, mm=None):
    """e00_{p,q+1}/e00_{p,q} = -i(2d-1)/(q+1), q = d-p."""
    dd = dd or d; mm = mm or m
    return -i*(2*dd-1)/(dd-mm+1)

# tau1_n, tau2_n: exact constants in Q(i,sqrt2)
def tau1(n):
    s = n+1
    return (-1)**s * ( (i-bq)**(-(s+1)) + (i+aq)**(-(s+1)) - i**(-(s+1)) )
def tau2(n):
    s = n+1
    return (-1)**s * ( (i-aq)**(-(s+1)) + (i+bq)**(-(s+1)) - i**(-(s+1)) )

# ---------------- zigzag top coefficients ----------------
def e_top(delta):
    """top Newton coefficient e(delta, floor(delta/2)) as a rational function of (d,m) times e00."""
    if delta % 2 == 0:
        s = delta//2
        return (-1)**s * A[s+1] / Pi(s)
    else:
        s = (delta-1)//2
        return (-1)**(s+1) * i * A[s+1] * (d-s-2) / Pi(s)

# ---------------- descending recursion ----------------
_e_cache = {}
def e_desc(delta, s):
    """candidate e(delta,s)/e00 for OFF-DIAGONAL keys (both parents off-diagonal),
    as a rational function of (d, m=p). Descends from the zigzag top via
    e(delta,s) = Lambda(delta,s) - e(delta,s+1)."""
    if s > delta//2: return sp.Integer(0)
    if s == delta//2: return e_top(delta)
    key = (delta, s)
    if key in _e_cache: return _e_cache[key]
    val = sp.expand(Lam(delta, s) - e_desc(delta, s+1))
    _e_cache[key] = val
    return val

def e_parent(delta, s, which):
    """e(delta,s)/e00 of the parent key, divided by e00 of the child:
    which='A': parent (p+1,q) of degree d+1: ratio includes rA0.
    which='B': parent (p,q+1): ratio includes rB0."""
    if which == 'A':
        return sp.expand(e_desc_at(delta, s, d+1, m+1) * rA0())
    else:
        return sp.expand(e_desc_at(delta, s, d+1, m) * rB0())

def e_desc_at(delta, s, dd, mm):
    """e(delta,s)/e00 evaluated at (d=dd, m=mm)."""
    if s > delta//2: return sp.Integer(0)
    if s == delta//2:
        if delta % 2 == 0:
            ss = delta//2
            out = (-1)**ss * A[ss+1] / Pi(ss, dd)
        else:
            ss = (delta-1)//2
            out = (-1)**(ss+1) * i * A[ss+1] * (dd-ss-2) / Pi(ss, dd)
        return sp.expand(out.subs({d: dd, m: mm}) if hasattr(out, 'subs') else out)
    # recursion at (dd,mm)
    return sp.expand(Lam_at(delta, s, dd, mm) - e_desc_at(delta, s+1, dd, mm))

def Lam_at(delta, s, dd, mm):
    """Lambda^(delta)(s) for key (m, d-m) = (mm, dd-mm), off-diagonal: parent feeds."""
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        cA = tau1(delta-dp) * (mm+1) / (2*dd-dp-1)
        cB = tau2(delta-dp) * (dd-mm+1) / (2*dd-dp-1)
        tot += cA * e_desc_at(dp, s, dd+1, mm+1) * rA0(dd, mm)
        tot += cB * e_desc_at(dp, s, dd+1, mm) * rB0(dd, mm)
    return sp.expand(tot)

def Lam(delta, s):
    return Lam_at(delta, s, d, m)
'''
open('tow2lib.py','w').write(lib)
print('tow2lib.py written')
