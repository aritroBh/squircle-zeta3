pipe = r'''"""cor6_pipe.py -- generic pipeline: true symbol of Omega(c0), delta_22, class."""
from fractions import Fraction as Fr
import json
from cor6_nf import NF, primes_above
from cor6_sunit import build_basis, exponent_vector
from cor6_sym import *
from cor6_ft import *
import sympy as sp

BASE = "<project>/research/"
D = json.load(open(BASE + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
W4TERMS = [(int(k.split(' ', 1)[0][2:]), eval(k.split(' ', 1)[1][1:]), Fr(v))
           for k, v in D if len(eval(k.split(' ', 1)[1][1:])) == 4]

class Corner:
    def __init__(self, label, minpoly, c0, S0=None, dimV=None, basis=None,
                 search_bound=3, search_denom=2, verbose=True):
        self.label = label
        self.F = NF(minpoly)
        F = self.F
        self.c0 = c0
        one = F.one()
        if S0 is None:
            S0 = set([2])
            for N in (F.norm(c0), F.norm(F.sub(one, c0)), F.norm(F.add(one, c0))):
                S0 |= set(sp.factorint(abs(N.numerator)).keys()) | set(sp.factorint(N.denominator).keys())
            S0 = sorted(p for p in S0 if p > 1)
        self.S0 = S0
        if dimV is None:
            r1 = F.d  # totally real
            nS = 0
            for p in S0:
                cnt, degs = primes_above(F.m, p)
                nS += len(degs)  # distinct primes
            dimV = r1 - 1 + nS
        self.dimV = dimV
        if basis is None:
            basis, pool = build_basis(F, S0, dimV, search_bound=search_bound,
                                      search_denom=search_denom, verbose=verbose)
        self.basis = basis
        self.letters_needed = {}
        self.zero, self.one = F.zero(), F.one()
        self.L = {0: F.zero(), 1: F.one(), -1: F.neg(F.one())}
    def vec(self, el):
        if el not in self.letters_needed:
            ev = exponent_vector(self.F, el, self.basis)
            if ev is None:
                raise ValueError(f"factor failed for {el} in {self.label}")
            self.letters_needed[el] = ev
        return self.letters_needed[el]
    def sym_H(self, word):
        wordF = tuple(self.L[a] for a in word)
        s = -1 if sum(1 for a in word if a == 1) % 2 else 1
        return sscale(s, sym_J(wordF, self.c0, self.vec, a0=self.zero,
                               sub=self.F.sub, iszero=lambda a: all(x == 0 for x in a)))
    def omega_symbol(self):
        SOm = {}
        for e, w, val in W4TERMS:
            ce = self.F.pow(self.c0, e)
            SOm = ft_add(SOm, ft_scale(self.F, self.F.scale(val, ce), self.sym_H(w)))
        return SOm
    def li4_symbol(self, x):
        vx, vomx = self.vec(x), self.vec(self.F.sub(self.one, x))
        return {(lneg(vomx), vx, vx, vx): Fr(1)}
'''
with open("cor6_pipe.py", "w") as f:
    f.write(pipe)

from fractions import Fraction as Fr
import importlib
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_pipe import Corner
from cor6_sym import W4Quotient
from cor6_ft import ft_coords, ft_delta22

# triangle corner: c = sqrt3/2, F = Q(sqrt3), theta = sqrt3, c0 = theta/2
C = Corner("c=sqrt3/2 (triangle, n=3)", [-3, 0, 1], (Fr(0), Fr(1, 2)))
print("basis:", C.basis)
for lbl, el in [("c", C.c0), ("1-c", C.F.sub(C.one, C.c0)), ("1+c", C.F.add(C.one, C.c0)), ("2", C.F.from_int(2))]:
    print(f"[{lbl}] = {C.vec(el)}")
SOm = C.omega_symbol()
print("words:", len(SOm))
d22 = ft_delta22(SOm, C.F.d, C.dimV)
print("delta_22(Omega(sqrt3/2)) nonzero components:", len(d22))
for k, v in sorted(d22.items()):
    print("   ", k, v)
