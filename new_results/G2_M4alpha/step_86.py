cwd = "<project>/research/08_theory/G2_M4alpha"
code = '''"""step_13.py -- TRUE number-level (3,2) pairings of the weight-5 part of B4 at corner points.

Adapts COR6 (08_theory/COR6_true_point_symbols) to weight 5. The weight-5 log-letter symbol of
M4(c0) is Omega_5(c0) = sum_{|w|=5} q c0^e [ S(H_w(c0)) + sgn S(H_{-w}(c0)) ]  (R4's weight-5 terms;
zeta parts and lower-weight terms have no weight-5 log-letter symbol). Letters: a certified S-unit
basis of the corner field (COR6 machinery; exponent vectors verified EXACTLY).

The (3,2) pairing follows w5_symbol.py's obstruction convention: with Q,T a basis of Lie_3(V)
(V 2-dimensional here: Q = [X,[X,Y]], T = [Y,[X,Y]]) and P = [X,Y] the Lambda^2 basis,
  obs_Q = sum_{u in supp Q, v in supp P} c_u d_v (s_{uv} - s_{vu}),   likewise obs_T.
Controls: (i) shuffle products are killed (exhaustive, all splits); (ii) at c0 = 1/sqrt(2) the
pairing must reproduce M4c's (84, -38 sqrt(2)) (w5_symbol.py, symbol.json) up to normalisation.
Then the same at c0 = 1/2 (hexagon) and c0 = 1/3. NAIVE vs TRUE: the naive pairing specialises the
function-level symbol (step_12) at the point; the TRUE one uses the S-unit letters directly.
"""
import sys, json, time
from fractions import Fraction as Fr
COR6 = "<project>/research/08_theory/COR6_true_point_symbols/"
sys.path.insert(0, COR6)
sys.path.insert(0, '.')
from cor6_nf import NF, primes_above
from cor6_sunit import build_basis, exponent_vector
from cor6_sym import sym_J, sscale, sadd
from cor6_ft import ft_add, ft_scale
import sympy as sp
import g2_fit as gf
from step_10 import load_terms
from step_12 import proj_symbol, S5_func, br, X as fX, Y as fY, Z as fZ


class W5Point:
    def __init__(self, label, minpoly, c0, verbose=False):
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
            nS = sum(len(primes_above(F.m, p)[1]) for p in self.S0)
            self.dimV = r1 - 1 + nS
            self.basis, _ = build_basis(F, self.S0, self.dimV, verbose=verbose)
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

    def omega5(self):
        S = {}
        for kind, w, e, q in load_terms(['hpl_fit_k4_(0,0,0).json']):
            if len(w) != 5:
                continue
            sgn = 1 if (e + gf.nz(w)) % 2 == 0 else -1
            ce = self.F.pow(self.c0, e)
            coef = self.F.scale(q, ce)
            S = ft_add(S, ft_scale(self.F, coef, self.sym_H(w)))
            if sgn == -1:
                coef = self.F.neg(coef)
            S = ft_add(S, ft_scale(self.F, coef, self.sym_H(gf.neg(w))))
        return S

    def letter_vectors(self):
        F, one = self.F, self.one
        return {lbl: self.vec(el) for lbl, el in
                [("c", self.c0), ("1-c", F.sub(one, self.c0)), ("1+c", F.add(one, self.c0))]}


def obstruction(S, dim):
    """[pairing([X,[X,Y]] x X^Y), pairing([Y,[X,Y]] x X^Y)] as NF elements."""
    X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
    Q = br(X, br(X, Y)); T = br(Y, br(X, Y)); P = br(X, Y)
    out = []
    for QQ in (Q, T):
        tot = None
        for u, cu in QQ.items():
            for v, dv in P.items():
                for word, sgn in ((u + v, 1), (v + u, -1)):
                    c = S.get(word)
                    if c is None:
                        continue
                    cc = tuple(sgn*cu*dv*x for x in c)
                    tot = cc if tot is None else tuple(a + b for a, b in zip(tot, cc))
        out.append(tot if tot is not None else (Fr(0),)*0)
    return out


def shuffle_control(P, n=30):
    import random, itertools
    random.seed(11)
    ok = True
    for t in range(n):
        lu = random.choice([1, 2, 3, 4]); lv = 5 - lu
        u = tuple(random.choice([-1, 0, 1]) for _ in range(lu))
        v = tuple(random.choice([-1, 0, 1]) for _ in range(lv))
        S = {}
        for w in gf_shuffles(u, v):
            for sw, cf in sadd(P.sym_H(w)).items():
                S[sw] = S.get(sw, Fr(0)) + cf
        # multiply nothing: pairing of the plain shuffle sum of symbols must vanish
        X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
        for QQ in (br(X, br(X, Y)), br(Y, br(X, Y))):
            PP = br(X, Y)
            tot = Fr(0)
            for uu, cu in QQ.items():
                for vv, dv in PP.items():
                    tot += cu*dv*(S.get(uu + vv, Fr(0)) - S.get(vv + uu, Fr(0)))
            if tot:
                ok = False
    return ok


def gf_shuffles(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + w for w in gf_shuffles(u[1:], v)] + [(v[0],) + w for w in gf_shuffles(u, v[1:])]


def naive_obstruction(pt):
    """specialise the FUNCTION-level weight-5 symbol (step_12) at the point's letter vectors."""
    lv = pt.letter_vectors()
    M = [lv["c"], lv["1-c"], lv["1+c"]]      # function letters [c],[1-c],[1+c] in the point basis
    S5 = S5_func()                            # {word over {1,2,3}: {e: Fr}}
    # evaluate Laurent polynomials at c0 (in F) and substitute letters
    dim = len(pt.basis)
    S = {}
    for word, d in S5.items():
        val = pt.F.zero()
        for e, cf in d.items():
            val = pt.F.add(val, pt.F.scale(cf, pt.F.pow(pt.c0, e)))
        newword = tuple(M[ix - 1] for ix in word)
        S[newword] = val
    # pair
    X = {(0,): Fr(1)}; Y = {(1,): Fr(1)}
    out = []
    for QQ in (br(X, br(X, Y)), br(Y, br(X, Y))):
        PP = br(X, Y)
        tot = None
        for u, cu in QQ.items():
            for v, dv in PP.items():
                for word, sgn in ((u + v, 1), (v + u, -1)):
                    c = S.get(word)
                    if c is None: continue
                    cc = tuple(sgn*cu*dv*x for x in c)
                    tot = cc if tot is None else tuple(a + b for a, b in zip(tot, cc))
        out.append(tot)
    return out


if __name__ == '__main__':
    t0 = time.time()
    # c0 = 1/sqrt(2) (square corner), F = Q(sqrt2)
    P4 = W5Point("c=1/sqrt2 (square)", [-2, 0, 1], (Fr(0), Fr(1, 2)))
    print(f"{P4.label}: field deg {P4.F.d}, S={P4.S0}, dimV={P4.dimV}, basis={P4.basis} ({time.time()-t0:.0f}s)", flush=True)
    S4 = P4.omega5()
    print(f"  weight-5 symbol: {len(S4)} tensor words ({time.time()-t0:.0f}s)", flush=True)
    obs4 = obstruction(S4, 2)
    print("  TRUE (3,2) pairing [Q, T]:", [[str(x) for x in c] for c in obs4], " (target [84, 0], [0, -38])")
    print("  letter vectors:", {k: [str(x) for x in v] for k, v in P4.letter_vectors().items()})
    n4 = naive_obstruction(P4)
    print("  NAIVE (specialised function) pairing:", [[str(x) for x in c] for c in n4])
    print("  shuffle control at this point:", shuffle_control(P4), flush=True)
    json.dump({'obs': [[str(x) for x in c] for c in obs4], 'naive': [[str(x) for x in c] for c in n4],
               'letters': {k: [str(x) for x in v] for k, v in P4.letter_vectors().items()}},
              open('point_obs_sqrt2.json', 'w'), indent=1)
    # c0 = 1/2 (hexagon corner, alpha = 2pi/3), F = Q
    P6 = W5Point("c=1/2 (hexagon)", [-1, 1], (Fr(1, 2),))
    print(f"{P6.label}: field deg {P6.F.d}, S={P6.S0}, dimV={P6.dimV}, basis={P6.basis} ({time.time()-t0:.0f}s)", flush=True)
    S6 = P6.omega5()
    obs6 = obstruction(S6, 2)
    print("  TRUE (3,2) pairing [Q, T]:", [[str(x) for x in c] for c in obs6])
    print("  letter vectors:", {k: [str(x) for x in v] for k, v in P6.letter_vectors().items()})
    n6 = naive_obstruction(P6)
    print("  NAIVE (specialised function) pairing:", [[str(x) for x in c] for c in n6])
    print("  shuffle control at this point:", shuffle_control(P6), flush=True)
    json.dump({'obs': [[str(x) for x in c] for c in obs6], 'naive': [[str(x) for x in c] for c in n6],
               'letters': {k: [str(x) for x in v] for k, v in P6.letter_vectors().items()}},
              open('point_obs_half.json', 'w'), indent=1)
    # c0 = 1/3, F = Q
    P3 = W5Point("c=1/3", [-1, 1], (Fr(1, 3),))
    S3 = P3.omega5()
    obs3 = obstruction(S3, 2)
    print(f"{P3.label}: S={P3.S0}, dimV={P3.dimV} ({time.time()-t0:.0f}s)")
    print("  TRUE (3,2) pairing [Q, T]:", [[str(x) for x in c] for c in obs3])
    n3 = naive_obstruction(P3)
    print("  NAIVE (specialised function) pairing:", [[str(x) for x in c] for c in n3])
    print("  shuffle control at this point:", shuffle_control(P3), flush=True)
    json.dump({'obs': [[str(x) for x in c] for c in obs3], 'naive': [[str(x) for x in c] for c in n3],
               'letters': {k: [str(x) for x in v] for k, v in P3.letter_vectors().items()}},
              open('point_obs_third.json', 'w'), indent=1)
'''
with open(cwd + "/step_13.py", "w") as f:
    f.write(code)
import ast; ast.parse(code)
print("written")
