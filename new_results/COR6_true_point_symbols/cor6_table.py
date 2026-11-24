"""cor6_table.py -- compute the full COR6 table: naive vs TRUE delta_22 at all points."""
from fractions import Fraction as Fr
import json, time
from cor6_nf import NF, primes_above
from cor6_sunit import build_basis, exponent_vector
from cor6_sym import *
from cor6_ft import *
import sympy as sp

BASE = "<project>/research/"
D = json.load(open(BASE + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
W4TERMS = [(int(k.split(' ', 1)[0][2:]), eval(k.split(' ', 1)[1][1:]), Fr(v))
           for k, v in D if len(eval(k.split(' ', 1)[1][1:])) == 4]

def vecQ(x):
    x = abs(Fr(x)); a = b = 0
    n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1, x
    return (Fr(a), Fr(b))

class Point:
    """A corner point c0 in a number field F (or Q with minpoly [-1,1])."""
    def __init__(self, label, minpoly, c0, basis="auto", verbose=False):
        self.label = label
        self.F = NF(minpoly)
        F = self.F
        self.c0 = c0
        one = F.one()
        S0 = set([2])
        for N in (F.norm(c0), F.norm(F.sub(one, c0)), F.norm(F.add(one, c0))):
            S0 |= set(sp.factorint(abs(N.numerator)).keys()) | set(sp.factorint(abs(N.denominator)).keys())
        self.S0 = sorted(p for p in S0 if p > 1)
        if F.d == 1:
            self.dimV = len(self.S0)
            self.basis = [F.from_int(p) for p in self.S0]
        else:
            r1 = F.d
            nS = 0
            for p in self.S0:
                cnt, degs = primes_above(F.m, p)
                nS += len(degs)
            self.dimV = r1 - 1 + nS
            if basis == "auto":
                self.basis, _ = build_basis(F, self.S0, self.dimV, verbose=verbose)
            elif basis == "letter-span":
                cands = [F.from_int(2), c0, F.sub(one, c0), F.add(one, c0)]
                bs = []
                for a in cands:
                    if not bs: bs.append(a); continue
                    if exponent_vector(F, a, bs, dps=120, S0=self.S0) is not None: continue
                    bs.append(a)
                self.basis = bs
            else:
                self.basis = basis
        self.letters = {}
        self.zero, self.one = F.zero(), F.one()
        self.L = {0: self.zero, 1: self.one, -1: F.neg(self.one)}
    def vec(self, el):
        if el not in self.letters:
            ev = exponent_vector(self.F, el, self.basis, S0=self.S0)
            if ev is None:
                raise ValueError(f"factor failed for {el} in {self.label}")
            self.letters[el] = ev
        return self.letters[el]
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
    def letter_vectors(self):
        F, one = self.F, self.one
        return {lbl: self.vec(el) for lbl, el in
                [("2", F.from_int(2)), ("c", self.c0), ("1-c", F.sub(one, self.c0)), ("1+c", F.add(one, self.c0))]}
    def letter_span_dim(self):
        lv = list(self.letter_vectors().values())
        rows = [list(map(Fr, v)) for v in lv]
        # exact rational rank
        if not rows: return 0
        M = rows
        R = 0
        for c in range(len(M[0])):
            p = next((r for r in range(R, len(M)) if M[r][c]), None)
            if p is None: continue
            M[R], M[p] = M[p], M[R]
            q = M[R][c]; M[R] = [x / q for x in M[R]]
            for r in range(len(M)):
                if r != R and M[r][c]:
                    M[r] = [a - b * M[r][c] for a, b in zip(M[r], M[R])]
            R += 1
        return R
    def naive_delta22(self):
        """COR4 function-level delta22 specialized: 
        (1-c^-1)([c]^[1-c])^([1-c]^[1+c]) + (-1-c^-1)([c]^[1+c])^([1-c]^[1+c])."""
        F = self.F
        vc = self.vec(self.c0)
        vm = self.vec(F.sub(self.one, self.c0))
        vp = self.vec(F.add(self.one, self.c0))
        dim = len(self.basis)
        out = {}
        def addw(u, v, x, y, coef):
            # (u^v)^(x^y) with NF coefficient coef
            w1 = wedge(u, v, dim); w2 = wedge(x, y, dim)
            for p1, c1 in w1.items():
                for p2, c2 in w2.items():
                    if p1 == p2: continue
                    key = (p1, p2) if p1 < p2 else (p2, p1)
                    sgn = Fr(1) if p1 < p2 else Fr(-1)
                    cc = F.scale(c1 * c2 * sgn, coef)
                    cur = out.get(key)
                    out[key] = cc if cur is None else tuple(a + b for a, b in zip(cur, cc))
        cinv = F.inv(self.c0)
        coef1 = F.sub(self.one, cinv)
        coef2 = F.neg(F.add(self.one, cinv))
        addw(vc, vm, vm, vp, coef1)
        addw(vc, vp, vm, vp, coef2)
        return {k: v for k, v in out.items() if any(v)}
    def true_delta22(self, SOm):
        return ft_delta22(SOm, self.F.d, len(self.basis))
