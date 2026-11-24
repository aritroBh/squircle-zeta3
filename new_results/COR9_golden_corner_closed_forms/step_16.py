import os, time
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
lib = '''"""cor9_lib.py -- shared machinery for COR9.
Reuses COR7 (cor7_lib: M3 weight split from G1's proved B3 term lists; quadrature with
precision boost; LLL relation finder) and cor7_pt/cor6 (exact S-unit + symbol machinery).
Adds: fundamental unit / norm-l elements via continued fractions (for large D census),
class number via reduced indefinite forms, Q(sqrt5) period bases.
"""
import sys
BASE = '<project>/research/'
sys.path.insert(0, BASE + '08_theory/COR7_classical_corners')
sys.path.insert(0, BASE + '08_theory/COR6_true_point_symbols')
import mpmath as mp
from fractions import Fraction as Fr
import cor7_lib as L7
from cor7_pt import QuadField, QuadPoint, factorint
from cor6_sym import W4Quotient, lneg
from cor6_ft import ft_coords, ft_delta22
from cor6_pipe import W4TERMS

Q2 = W4Quotient(2)

def li4_coords(P, x):
    vx, vo = P.vec(x), P.vec(P.F.sub(P.F.one(), x))
    return Q2.coords({(lneg(vo), vx, vx, vx): Fr(1)})

def omega_class(P):
    return ft_coords(Q2, P.omega_symbol(W4TERMS), 2)

def rankQ(vecs):
    M = [list(v) for v in vecs]; R = 0
    for c in range(len(M[0])):
        p = next((r for r in range(R, len(M)) if M[r][c]), None)
        if p is None: continue
        M[R], M[p] = M[p], M[R]
        q = M[R][c]; M[R] = [x/q for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][c]:
                M[r] = [a - b*M[R][c] for a, b in zip(M[r], M[R])]
        R += 1
    return R

# ---------- exact fundamental unit and norm-l elements via continued fractions ----------
def cf_sqrt(D):
    """continued fraction of sqrt(D): (a0, period)."""
    a0 = int(D**0.5)
    if a0*a0 == D: return a0, []
    m, d, a = 0, 1, a0
    per = []
    while True:
        m = d*a - m
        d = (D - m*m)//d
        a = (a0 + m)//d
        per.append(a)
        if a == 2*a0: break
    return a0, per

def pell_min(D, N_target):
    """minimal solution x^2 - D y^2 = N_target (N_target in {+1,-1,+4,-4}) via CF of sqrt D.
    Returns (x, y) positive or None."""
    a0, per = cf_sqrt(D)
    convs = []
    h1, h2, k1, k2 = 1, 0, 0, 1
    seq = [a0] + per*4
    for a in seq[:4*len(per)+2]:
        h1, h2 = a*h1 + h2, h1
        k1, k2 = a*k1 + k2, k1
        if h2*h2 - D*k2*k2 == N_target:
            return (h2, k2)
    return None

def fund_unit(D):
    """fundamental unit of Q(sqrt D) as (Fr(a,e), Fr(b,e)) = (a+b sqrt D)/e, e in {1,2}, > 1."""
    if D % 4 == 1:
        s = pell_min(D, -4) or pell_min(D, 4)
        if s: return (Fr(s[0], 2), Fr(s[1], 2))
    s = pell_min(D, -1) or pell_min(D, 1)
    return (Fr(s[0]), Fr(s[1]))

def norm_l_element(D, l):
    """element of O_F with |N| = l, as (a,b,e), via CF reduction on x^2 - D y^2 = +-l, +-4l;
    returns (Fr(a,e), Fr(b,e)) or None. Uses sympy diop_DN (Lagrange)."""
    from sympy.solvers.diophantine.diophantine import diop_DN
    best = None
    for tgt in ([-l, l] + ([-4*l, 4*l] if D % 4 == 1 else [])):
        try:
            sols = diop_DN(D, tgt)
        except Exception:
            sols = []
        e = 2 if abs(tgt) == 4*l else 1
        for x, y in sols:
            if e == 2 and (x % 2 != y % 2):
                continue
            el = (Fr(x, e), Fr(y, e))
            # verify norm exactly
            nn = el[0]*el[0] - D*el[1]*el[1]
            if abs(nn) != l: continue
            v = abs(float(el[0]) + float(el[1])*D**0.5)
            if v > 0:
                lg = abs(mp.log(mp.mpf(v)))
                if best is None or lg < best[0]:
                    best = (lg, el)
    return best[1] if best else None

# ---------- class number of real quadratic field via reduced indefinite forms ----------
def disc_F(D):
    return D if D % 4 == 1 else 4*D

def h_plus(disc):
    """narrow class number of the order of discriminant disc > 0 (non-square):
    number of cycles of reduced primitive indefinite forms. Cohen CCANT 5.6/5.7."""
    import math
    forms = []
    b0 = int(math.isqrt(disc))
    # enumerate reduced primitive forms (a,b,c), b^2 - 4ac = disc
    red = []
    for b in range(1, b0+1):
        if (disc - b*b) % 4: continue
        for a in range(1, b0+1):
            if (b*b - disc) % (4*a): continue
            c = (b*b - disc)//(4*a)
            if c == 0: continue
            g = math.gcd(math.gcd(a, b), abs(c))
            if g != 1: continue
            # reduced: |sqrt(disc)-2|a|| < b < sqrt(disc)
            if not (abs(math.isqrt(disc) - 2*abs(a)) < b ): continue
            red.append((a, b, c))
    # rho operator
    import math as m2
    sq = m2.isqrt(disc)
    def rho(f):
        a, b, c = f
        # xi = floor((b + sq) / (2c))
        xi = (b + sq)//(2*abs(c)) 
        # adjust sign convention: use rho(a,b,c) = (c, -b + 2 c t, a - b t + c t^2), t = floor((b+sq)/(2c)) with c signed
        t = (b + sq)//(2*c) if c > 0 else -((-b + sq)//(-2*c)) if False else None
        # do floor division properly
        t = (b + sq)//(2*c)  # python floor division works for negative c
        return (c, -b + 2*c*t, a - b*t + c*t*t)
    # build cycles
    seen = {}
    ncyc = 0
    for f in red:
        if f in seen: continue
        ncyc += 1
        g = f
        while g not in seen:
            seen[g] = ncyc
            g = rho(g)
    return ncyc

def class_number(D):
    """ordinary class number of Q(sqrt D), D squarefree > 1."""
    disc = disc_F(D)
    hp = h_plus(disc)
    # narrow = ordinary iff N(eps) = -1
    eps = fund_unit(D)
    N = eps[0]*eps[0] - D*eps[1]*eps[1]
    return hp if N == -1 else hp//2
'''
open(WD+'cor9_lib.py','w').write(lib)
print("lib written")
# quick test of fund_unit and class numbers
sys.path.insert(0, WD)
from cor9_lib import fund_unit, class_number, pell_min
for D in [2,3,5,6,7,14,46,61,94]:
    print(D, fund_unit(D), class_number(D))
