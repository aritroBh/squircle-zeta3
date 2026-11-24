lib = '''"""cor7_pt.py -- exact letter-vector machinery for quadratic corners with small S,
using searched-and-exactly-certified generators (no log-lattice inversion).
For a real quadratic field F = Q(sqrt D) and S = {2}:
  2 ramified: S-units = <-1, p, eps> (p^2 = 2*unit)
  2 inert:    S-units = <-1, 2, eps>
  2 split:    S-units = <-1, 2, p, eps> (p' = -2/p up to units)
vec(x) = exponent tuple over the generators, found by bounded exact search;
all searches are certified by exact NF multiplication.
"""
from fractions import Fraction as Fr
import itertools, sys
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from cor6_nf import NF
from cor6_sym import sym_J, sscale, delta22, sadd

def normD(D, a, b, e=1):
    return (a*a - D*b*b) / (e*e) if e != 1 else a*a - D*b*b

def find_unit(F, D, bound=80):
    """fundamental unit: N = +-1, smallest |log|>0 (search plain and half-integers)."""
    best = None
    cands = []
    for e in (1, 2):
        if e == 2 and D % 4 != 1:
            continue
        for a in range(-bound, bound+1):
            for b in range(1, bound+1):
                if (a*a - D*b*b) * (1 if e == 1 else 1) != 0:
                    if abs(a*a - D*b*b) == e*e:
                        cands.append((Fr(a, e), Fr(b, e)))
    import mpmath as mp
    mp.mp.dps = 60
    rD = mp.sqrt(D)
    for cf in cands:
        v = abs(float(cf[0]) + float(cf[1]) * float(rD))
        if v > 0:
            lg = abs(mp.log(v))
            if lg > 1e-9 and (best is None or lg < best[0]):
                best = (lg, cf)
    return best[1] if best else None

def find_p(F, D, bound=80):
    """element with N = +-2 (then p generates a prime above 2), smallest log."""
    best = None
    import mpmath as mp
    mp.mp.dps = 60
    rD = mp.sqrt(D)
    for e in (1, 2):
        if e == 2 and D % 4 != 1:
            continue
        for a in range(-bound, bound+1):
            for b in range(1, bound+1):
                if abs(a*a - D*b*b) == 2*e*e:
                    v = abs(a + b*float(rD))/e
                    lg = abs(mp.log(v)) if v > 0 else mp.inf
                    if best is None or lg < best[0]:
                        best = (lg, (Fr(a, e), Fr(b, e)))
    return best[1] if best else None

class QuadPoint:
    """corner point c0 in Q(sqrt D), letters in the {2}-unit group."""
    def __init__(self, label, D, c0, K=12):
        self.label = label
        self.D = D
        self.F = NF([-D, 0, 1])
        self.c0 = c0
        F = self.F
        self.eps = find_unit(F, D)
        self.p = find_p(F, D)   # None if 2 inert
        # generator list: (2 or p-ramified or 2-inert), eps, [p]
        ramified = (D % 8 in (2, 3)) or (D % 4 == 0)
        split = (D % 8 == 1)
        inert = not ramified and not split
        gens = []
        if inert:
            gens = [F.from_int(2), self.eps]
        elif ramified:
            gens = [self.p, self.eps]
            # verify p^2 = 2*unit: N(p)^2/N(2): |N(p)| = 2: p^2/2 = unit
        else:  # split
            gens = [F.from_int(2), self.p, self.eps]
        self.gens = gens
        self.split, self.ramified, self.inert = split, ramified, inert
        self.cache = {F.one(): (0,)*len(gens), F.neg(F.one()): (0,)*len(gens)}
        self.K = K
        zero = F.zero(); one = F.one()
        self.L = {0: zero, 1: one, -1: F.neg(one)}
    def vec_search(self, el):
        F = self.F
        if el in self.cache:
            return self.cache[el]
        K = self.K
        n = len(self.gens)
        for tup in itertools.product(range(-K, K+1), repeat=n):
            v = F.one()
            for e_, g in zip(tup, self.gens):
                if e_:
                    v = F.mul(v, F.pow(g, e_))
            if F.eq(v, el):
                self.cache[el] = tup
                return tup
        return None
    def vec(self, el):
        v = self.vec_search(el)
        if v is not None:
            return v
        # try multiplicative build-up from known cache (cheap closure)
        for e1, v1 in list(self.cache.items()):
            q = self.F.mul(el, self.F.inv(e1))
            vq = self.vec_search(q) if q not in self.cache else self.cache[q]
            if vq is not None:
                out = tuple(a - b for a, b in zip(vq, v1))
                self.cache[el] = out
                return out
        raise ValueError(f"vec failed for {el} in {self.label}")
    def sym_H(self, word):
        wordF = tuple(self.L[a] for a in word)
        s = -1 if sum(1 for a in word if a == 1) % 2 else 1
        return sscale(s, sym_J(wordF, self.c0, self.vec, a0=self.F.zero(),
                               sub=self.F.sub, iszero=lambda a: all(x == 0 for x in a)))
    def letter_span_rank(self):
        F, one = self.F, self.F.one()
        vecs = [self.vec(F.from_int(2)), self.vec(self.c0),
                self.vec(F.sub(one, self.c0)), self.vec(F.add(one, self.c0))]
        M = [list(map(Fr, v)) for v in vecs]
        R = 0
        for c in range(len(M[0])):
            p_ = next((r for r in range(R, len(M)) if M[r][c]), None)
            if p_ is None: continue
            M[R], M[p_] = M[p_], M[R]
            q = M[R][c]; M[R] = [x/q for x in M[R]]
            for r in range(len(M)):
                if r != R and M[r][c]:
                    M[r] = [a - b*M[R][c] for a, b in zip(M[r], M[R])]
            R += 1
        return R, vecs

def omega_symbol_pt(P, W4TERMS):
    """symbol of Omega(c0) over the letter vectors (Fr tensors), coefficients c0^e exact."""
    import sys
    sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
    from cor6_ft import ft_add, ft_scale
    SOm = {}
    for e, w, val in W4TERMS:
        ce = P.F.pow(P.c0, e)
        SOm = ft_add(SOm, ft_scale(P.F, P.F.scale(val, ce), P.sym_H(w)))
    return SOm
'''
open('cor7_pt.py', 'w').write(lib)
print("cor7_pt.py written")
