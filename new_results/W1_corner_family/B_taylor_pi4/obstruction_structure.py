"""obstruction_structure.py -- EXACT (symbol level).
(1) Write delta_{2,2}(A_e), e = -1, 0, as wedges of Bloch-group images rho{z}_2 := z ^ (1-z) of natural dilogs.
(2) Pull back to tau = tan(alpha/4) (c = (1-tau^2)/(1+tau^2), s = sin(alpha/2) = 2tau/(1+tau^2));
    the involution iota: tau -> (1-tau)/(1+tau) (alpha -> pi - alpha, c <-> s) fixes tau = sqrt2 - 1 (the square
    corner, c = s = 1/sqrt2).  Compute the rank of the obstructions of A_{-1}(c), A_0(c), A_{-1}(s), A_0(s), and of
    the classical A_1(c), A_1(s), in Lambda^2(Lambda^2 <2, tau, 1-tau, 1+tau, 1+tau^2>).
    Rank 4 means: no linear combination (with any rational-function prefactors) of these four non-classical
    components is classical, so no prefactor reshuffling at the fixed point c = s produces a classical family.
(3) Specialise the obstruction to c = 1/sqrt2: the letters 2, c, 1-c, 1+c map into <sqrt2, a> (a = 1+sqrt2), a rank-2
    group, where Lambda^2(Lambda^2) = 0.  Print the specialised delta_{2,2} (must be 0) and the specialised
    Lambda^2 images of the two B_2 factors (they become proportional).
"""
import itertools, sys
from fractions import Fraction as Fr
import sympy as sp
sys.path.insert(0, '.')
import hpl
from b2b2_test import load, lin, delta22, PAIRS

R3 = load('hpl_fit_k3_(0,0,0).json')
comp = {}
for e in (-1, 0, 1):
    S = {}
    for ee, w, v in R3:
        if ee == e and len(w) == 4:
            S = lin((1, S), (v, hpl.symbol(w)))
    comp[e] = S

# ---- (1) B2 images in Lambda^2 <2, x, 1-x, 1+x>  (vectors over PAIRS)
def wedge_vec(u, v):
    """u, v exponent vectors over the letter basis; returns vector over PAIRS of u^v."""
    out = [Fr(0)]*len(PAIRS)
    for i in range(len(u)):
        for j in range(len(v)):
            if i == j or u[i] == 0 or v[j] == 0:
                continue
            if i < j:
                out[PAIRS.index((i, j))] += u[i]*v[j]
            else:
                out[PAIRS.index((j, i))] -= u[i]*v[j]
    return out

# letters: index 0 = 2, 1 = x, 2 = 1-x, 3 = 1+x  (mod torsion)
Z = {  # z : (vector of z, vector of 1-z)
    'x': ((0, 1, 0, 0), (0, 0, 1, 0)),
    '-x': ((0, 1, 0, 0), (0, 0, 0, 1)),
    'x^2': ((0, 2, 0, 0), (0, 0, 1, 1)),
    '(1-x)/2': ((-1, 0, 1, 0), (-1, 0, 0, 1)),
    '-(1-x)/(1+x)': ((0, 0, 1, -1), (1, 0, 0, -1)),
    '(1-x)/(1+x)': ((0, 0, 1, -1), (1, 1, 0, -1)),
}
rho = {k: wedge_vec(z, omz) for k, (z, omz) in Z.items()}


def outer_wedge(X, Y):
    out = {}
    for P in range(len(PAIRS)):
        for Q in range(len(PAIRS)):
            if X[P] == 0 or Y[Q] == 0 or P == Q:
                continue
            cf = X[P]*Y[Q]
            if P < Q:
                out[(P, Q)] = out.get((P, Q), 0) + cf
            else:
                out[(Q, P)] = out.get((Q, P), 0) - cf
    return {k: v for k, v in out.items() if v != 0}

Y = rho['-(1-x)/(1+x)']
Qa = rho['x^2']
Qb = [p - q for p, q in zip(rho['x'], rho['-x'])]
for e in (-1, 0):
    D = delta22(comp[e])
    for name, Qv in (('rho{x^2}', Qa), ('rho{x} - rho{-x}', Qb)):
        W = outer_wedge(Y, Qv)
        # find lambda with D = lambda * W
        ratios = {D.get(k, 0)/W[k] for k in W} if W else set()
        if W and len(ratios) == 1 and set(D) == set(W):
            print(f'(1) delta22(A_{e}) = {ratios.pop()} * rho{{-(1-c)/(1+c)}} ^ {name}')

# ---- (2) pullback to tau
t = sp.Symbol('t')
TL = [sp.Integer(2), t, 1 - t, 1 + t, 1 + t**2]   # tau-letters


def factor_vec(expr):
    expr = sp.factor(sp.cancel(expr))
    vec = [0]*5
    num, den = sp.fraction(expr)
    for poly, sg in ((num, 1), (den, -1)):
        c, facs = sp.factor_list(poly, t)
        # constant: powers of 2 only (mod sign)
        cc = abs(sp.Rational(c))
        k2 = 0
        while cc % 2 == 0 and cc != 0:
            cc /= 2; k2 += 1
        assert cc == 1 or (1/cc) % 1 == 0, expr
        while cc != 1 and (1/cc) % 2 == 0:
            cc *= 2; k2 -= 1
        assert cc == 1, (expr, c)
        vec[0] += sg*k2
        for f, p in facs:
            f = sp.expand(f)
            idx = None
            for j, L in enumerate(TL[1:], 1):
                if sp.expand(f - L) == 0 or sp.expand(f + L) == 0:
                    idx = j
            assert idx is not None, (expr, f)
            vec[idx] += sg*p
    return vec

cfun = (1 - t**2)/(1 + t**2)
sfun = 2*t/(1 + t**2)
def letter_map(xf):
    return [factor_vec(sp.Integer(2)), factor_vec(xf), factor_vec(1 - xf), factor_vec(1 + xf)]

def pullback(S, xf):
    lm = letter_map(xf)
    out = {}
    for word, cf in S.items():
        vecs = [lm[i] for i in word]
        for idx in itertools.product(*[[(k, e) for k, e in enumerate(v) if e] for v in vecs]):
            key = tuple(k for k, _ in idx)
            val = cf
            for _, e in idx:
                val *= e
            out[key] = out.get(key, 0) + val
    return {k: v for k, v in out.items() if v != 0}

PAIRS5 = [(i, j) for i in range(5) for j in range(i + 1, 5)]
def delta22_5(S):
    out = {}
    for (a, b, c, d), cf in S.items():
        if a == b or c == d:
            continue
        s1 = 1 if a < b else -1; P = PAIRS5.index((min(a, b), max(a, b)))
        s2 = 1 if c < d else -1; Q = PAIRS5.index((min(c, d), max(c, d)))
        if P == Q:
            continue
        sg = s1*s2
        if P > Q:
            P, Q, sg = Q, P, -sg
        out[(P, Q)] = out.get((P, Q), 0) + sg*cf
    return {k: v for k, v in out.items() if v != 0}

obs = {}
for e in (-1, 0, 1):
    for nm, xf in (('c', cfun), ('s', sfun)):
        obs[(e, nm)] = delta22_5(pullback(comp[e], xf))
        print(f'(2) delta22 of A_{e}({nm}(tau)) in the tau alphabet: {len(obs[(e, nm)])} nonzero components')
keys = [(-1, 'c'), (0, 'c'), (-1, 's'), (0, 's')]
rows = sorted(set().union(*[set(obs[k]) for k in keys]))
Mx = sp.Matrix([[obs[k].get(r, 0) for k in keys] for r in rows])
print('(2) rank of the obstructions of A_{-1}(c), A_0(c), A_{-1}(s), A_0(s):', Mx.rank())

# ---- (3) specialise at c = 1/sqrt2: letters 2, c, 1-c, 1+c -> exponent vectors over (sqrt2, a)
spec = {0: (2, 0), 1: (-1, 0), 2: (-1, -1), 3: (-1, 1)}   # 2 = sqrt2^2, c = sqrt2^-1, 1-c = b/sqrt2, 1+c = a/sqrt2
def spec_wedge(vecpairs):
    # Lambda^2 <sqrt2, a> is 1-dim: (p1,q1)^(p2,q2) = (p1 q2 - q1 p2) e,  e = sqrt2 ^ a
    return vecpairs
def spec_rho(z_vec, omz_vec):
    p = [sum(z_vec[i]*spec[i][0] for i in range(4)), sum(z_vec[i]*spec[i][1] for i in range(4))]
    q = [sum(omz_vec[i]*spec[i][0] for i in range(4)), sum(omz_vec[i]*spec[i][1] for i in range(4))]
    return p[0]*q[1] - p[1]*q[0]
print('(3) at c = 1/sqrt2, Lambda^2-images (multiples of e = sqrt2 ^ a):')
for k, (z, omz) in Z.items():
    print(f'      rho{{{k}}} -> {spec_rho(z, omz)} e')
print('    so rho{-(1-c)/(1+c)} ^ rho{c^2} -> 0 and rho{-(1-c)/(1+c)} ^ (rho{c}-rho{-c}) -> 0 (Lambda^2 of a line).')
