"""OCT9d common: exact symbol targets, candidate pools, and evaluation helpers."""
import os, sys, pickle
base = "<project>/research/08_theory/"
O9B = base+"OCT9b_octagon_debug"
for p in (O9B, base+"OCT9c_octagon_close"):
    if p not in sys.path: sys.path.insert(0, p)
import sympy as sp
from fractions import Fraction as Fr
from syminfra import S_Li3, S_mLi2logp, S_logprod, SymbolSpace, symbol_of_expr
from oct9c_common import load_all, get_engine, rhs_terms_235, rhs_terms_2357

SS235 = SymbolSpace((2,3,5))
SS2357 = SymbolSpace((2,3,5,7))
D = load_all()
Gargs, Dargs = rhs_terms_235()
sols7, negw7, logprods7 = rhs_terms_2357()

HERE = os.path.dirname(os.path.abspath(__file__))
def TW3():
    """target symbol vectors (list of sympy Rationals), e=0..3"""
    return pickle.load(open(os.path.join(HERE, "TW3.pkl"), "rb"))

# ---------- candidate pool construction ----------
def factor_ok(q, primes):
    """is |q| a product of the given primes?"""
    q = abs(Fr(q))
    n, d = q.numerator, q.denominator
    for p in primes:
        while n % p == 0 and n > 1: n //= p
        while d % p == 0 and d > 1: d //= p
    return n == 1 and d == 1

def orbit(z):
    z = Fr(z)
    out = set()
    pts = [z, 1-z]
    try: pts.append(1/z)
    except ZeroDivisionError: pass
    for w in list(pts):
        out.add(w)
    # full S3 orbit
    s = {z}
    if z != 0: s.add(1/z)
    s.add(1-z)
    if z != 1 and z != 0:
        s.add(1/(1-z)); s.add(z/(z-1)); s.add((z-1)/z)
    return s

def good(z, primes):
    """z and 1-z both S-units (so S(Li3(z)) lives in the alphabet); z != 0,1"""
    z = Fr(z)
    if z in (0, 1): return False
    return factor_ok(z, primes) and factor_ok(1-z, primes)

# edge points: x1=-k1, x2=-k2, b1=k1/(1+k1), b2=k2/(1+k2)
K12 = {0: (Fr(3,2), Fr(-7,12)), 1: (Fr(1,3), Fr(1,9)),
       2: (Fr(5,3), Fr(-5,9)), 3: (Fr(1,4), Fr(1,24))}
def edge_points(e):
    k1, k2 = K12[e]
    return [ -k1, -k2, k1/(1+k1), k2/(1+k2) ]

SMALL235 = [Fr(1,2), Fr(1,3), Fr(1,4), Fr(2,3), Fr(3,4), Fr(1,5), Fr(2,5), Fr(4,5),
            Fr(1,6), Fr(5,6), Fr(1,9), Fr(1,10), Fr(9,10), Fr(5,9), Fr(4,9), Fr(8,9),
            Fr(1,25), Fr(24,25), Fr(1,16), Fr(15,16), Fr(1,81), Fr(80,81)]
SMALL7 = [Fr(1,7), Fr(6,7), Fr(1,8), Fr(7,8), Fr(1,49), Fr(48,49), Fr(7,15), Fr(7,12),
          Fr(5,12), Fr(1,50), Fr(49,50)]

def pool(edges, primes, small, add_edge_orbits=True):
    """candidate Li3 arguments: edge points + S3 orbits + small units + their orbits."""
    P = set()
    for e in edges:
        for z in edge_points(e):
            P |= orbit(z)
    for z in small:
        P |= orbit(z)
    P = {z for z in P if good(z, primes)}
    return sorted(P, key=lambda z: (abs(z).denominator, abs(z).numerator))

def col_Li3(g, SS): return SS.coord(S_Li3(g))
def col_prod(w, p, SS): return SS.coord(S_mLi2logp(w, p))
def col_logprod(trip, SS): return SS.coord(S_logprod(trip))
