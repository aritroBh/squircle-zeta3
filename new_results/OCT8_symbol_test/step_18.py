import os
code = r'''
"""OCT8: symbol machinery for Li_{2,1} reducibility tests (Goncharov criterion).
Symbols live in V^{otimes 3}, V = F^x (x) Q. For F=Q: letters = primes.
S(Li_3(z)) = -(1-z) (x) z (x) z.
S(Li_{2,1}(x,y)) = (xy-1)(x)(xy)(x)(y-1) + (xy-1)(x)[x(y-1)/(x-1)](x)x + (y-1)(x)(x-1)(x)x
   (verified: integrable; stuffle with Li_{1,2}(y,x)+Li_3(xy); Landen for Li_3).
Decomp = span of shuffles  a sh (b(x)c) = a(x)b(x)c + b(x)a(x)c + b(x)c(x)a .
Reducible-to-Li3 (mod products & constants)  iff  Pi S(target) in span{Pi S(Li_3(z))} in Q = V3/Decomp.
"""
from fractions import Fraction as Fr
import sympy as sp
import itertools

# ---------- vectors (elements of V) ----------
def vadd(*vs):
    r = {}
    for v in vs:
        for k, c in v.items():
            r[k] = r.get(k, Fr(0)) + c
            if r[k] == 0: del r[k]
    return r
def vmul(v, a):
    a = Fr(a); return {k: c*a for k, c in v.items() if c*a != 0}

def vec_of_rational(q):
    q = Fr(q); n, d = abs(q.numerator), q.denominator
    out = {}
    for m, sgn in ((n,1),(d,-1)):
        for p, e in sp.factorint(m).items():
            out[p] = out.get(p, Fr(0)) + sgn*e
    return out

# ---------- symbols ----------
def tensor3(u, v, w):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            for c, cc in w.items():
                S[(a,b,c)] = S.get((a,b,c), Fr(0)) + ca*cb*cc
    return {k: x for k, x in S.items() if x != 0}
def sadd(*Ss):
    r = {}
    for S in Ss:
        for k, c in S.items():
            r[k] = r.get(k, Fr(0)) + c
            if r[k] == 0: del r[k]
    return r
def smul(S, a):
    a = Fr(a); return {k: c*a for k, c in S.items() if c*a != 0}

def S_Li3(z, vecfun):
    """symbol of Li_3(z), z algebraic; vecfun: element -> vec"""
    zv = vecfun(z); omzv = vecfun(1-z)
    return tensor3(vmul(omzv, -1), zv, zv)

def S_Li21(x, y, vecfun):
    L1 = vecfun(x*y - 1)
    L2 = vecfun(x*y)
    L3 = vecfun(y - 1)
    L4 = vecfun(x*(y-1)/(x-1))
    L5 = vecfun(x - 1)
    L6 = vecfun(x)
    return sadd(tensor3(L1, L2, L3), tensor3(L1, L4, L6), tensor3(L3, L5, L6))

# ---------- quotient Q = V^3 / Decomp (exact rational linear algebra) ----------
class Quotient:
    def __init__(self, letters):
        self.letters = list(letters)
        self.idx = {l:i for i,l in enumerate(self.letters)}
        self.dim3 = len(self.letters)**3
        self.basis3 = list(itertools.product(self.letters, repeat=3))
        self.bidx = {t:i for i,t in enumerate(self.basis3)}
        # decomposables
        D = []
        for a in self.letters:
            for b in self.letters:
                for c in self.letters:
                    v = sp.zeros(1, self.dim3)
                    for t in [(a,b,c),(b,a,c),(b,c,a)]:
                        v[self.bidx[t]] += 1
                    D.append(v)
        self.D = D
        DM = sp.Matrix.hstack(*[d.T for d in D])
        self.rank_D = DM.rank()
    def coord(self, S):
        v = sp.zeros(1, self.dim3)
        for t, c in S.items():
            v[self.bidx[t]] += sp.Rational(c.numerator, c.denominator)
        return v
    def in_span(self, target, gens):
        """is target in span(gens)+Decomp ? returns (bool, coeffs or None)"""
        M = sp.Matrix.hstack(*[g.T for g in gens] + [d.T for d in self.D])
        r0 = M.rank()
        M2 = sp.Matrix.hstack(M, target.T)
        if M2.rank() != r0:
            return False, None
        # solve M c = target
        c = M.gauss_jordan_solve(target.T)[0]
        return True, [sp.Rational(x) for x in c[:len(gens)]]
    def witness(self, target, gens):
        """row w: w.D=0, w.gens=0, w.target != 0 (certificate of non-membership)."""
        M = sp.Matrix.hstack(*[g.T for g in gens] + [d.T for d in self.D])  # dim3 x m
        NS = M.T.nullspace()  # rows orthogonal to columns of M
        for w in NS:
            if w*target.T != sp.zeros(1,1):
                return w, (w*target.T)[0]
        return None, None
'''
open("oct8.py","w").write(code)
print("oct8.py written", len(code))
