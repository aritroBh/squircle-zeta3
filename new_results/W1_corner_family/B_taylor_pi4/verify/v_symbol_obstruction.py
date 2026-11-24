"""v_symbol_obstruction.py -- verifier's independent check of the function-level classicality obstruction.

Own code (does not import the solver's hpl.py / b2b2_test.py).
Function-level symbol of an HPL, from its differential d H_{a,w}(c) = f_a(c) H_w(c) dc:
    S(H_{a_1..a_k}) = l_{a_k} (x) ... (x) l_{a_1},  l_0 = [c], l_1 = -[1-c], l_{-1} = [1+c]
(letters modulo constants; basis e0 = [c], e1 = [1-c], e2 = [1+c]).
Obstruction operator  O(a(x)b(x)c(x)d) = (a^b)^(c^d)  in Lambda^2(Lambda^2 V).
It annihilates the symbol of Li_4(f) for every f (last two slots equal) and every shuffle product of lower
weights (checked below on random symbols), so O(S) != 0 excludes any classical expression (Li_4's and
products, C-coefficients) for the length-4 part.
Checks:
 (1) O(A_{-1}), O(A_0), O(A_1) with A_e the weight-4 HPL part of G1's R_3 with prefactor c^e;
 (2) comparison with the solver's rho-forms:  O(A_-1) = -1/2 rho{-(1-c)/(1+c)} ^ rho{c^2},
     O(A_0) = rho{-(1-c)/(1+c)} ^ (rho{c} - rho{-c}),  rho{z} = [z]^[1-z]  (function level);
 (3) pullback to tau (c = (1-t^2)/(1+t^2), s = 2t/(1+t^2)): rank of O of A_-1(c), A_0(c), A_-1(s), A_0(s);
 (4) specialisation of the letters at c0 = 1/sqrt2 into <sqrt2, a> (constant-retaining letters incl. 2);
 (5) independence of the two obstruction vectors O(A_-1), O(A_0) (so no combination r A_-1 + s A_0 with
     r, s not both 0 is classical).
"""
import itertools, json, re, random
from fractions import Fraction as Fr
import sympy as sp

G1 = 'research/10_frontier/G1_zeta3_any_angle/'
ELL = {0: {0: 1}, 1: {1: -1}, -1: {2: 1}}     # letter vectors over (c, 1-c, 1+c)


def load(fname):
    out = []
    for nm, v in json.load(open(G1 + fname)):
        m = re.fullmatch(r'c\^(-?\d+) H\[([-0-9, ]+)\]', nm)
        if m:
            out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), Fr(v)))
    return out


def tensor(vecs, cf=Fr(1)):
    """expand a tensor of letter vectors into {tuple of basis indices: coefficient}"""
    out = {}
    for combo in itertools.product(*[list(v.items()) for v in vecs]):
        key = tuple(k for k, _ in combo)
        val = cf
        for _, e in combo:
            val *= e
        out[key] = out.get(key, 0) + val
    return {k: v for k, v in out.items() if v != 0}


def hsym(w):
    return tensor([ELL[a] for a in reversed(w)])


def addS(S, T, cf=1):
    for k, v in T.items():
        S[k] = S.get(k, 0) + cf*v
        if S[k] == 0:
            del S[k]
    return S


def wedge2(i, j):
    """basis wedge e_i ^ e_j -> (sign, (min,max)) or None"""
    if i == j:
        return None
    return (1, (i, j)) if i < j else (-1, (j, i))


def O(S):
    out = {}
    for (a, b, c, d), cf in S.items():
        x, y = wedge2(a, b), wedge2(c, d)
        if x is None or y is None or x[1] == y[1]:
            continue
        sgn = x[0]*y[0]
        P, Q = x[1], y[1]
        if P > Q:
            P, Q, sgn = Q, P, -sgn
        out[(P, Q)] = out.get((P, Q), 0) + sgn*cf
        if out[(P, Q)] == 0:
            del out[(P, Q)]
    return out


def shuffle(S1, S2):
    out = {}
    for w1, c1 in S1.items():
        for w2, c2 in S2.items():
            n1, n2 = len(w1), len(w2)
            for pos in itertools.combinations(range(n1 + n2), n1):
                it1, it2 = iter(w1), iter(w2)
                word = tuple(next(it1) if k in pos else next(it2) for k in range(n1 + n2))
                out[word] = out.get(word, 0) + c1*c2
    return {k: v for k, v in out.items() if v != 0}


def vec_wedge(u, v):
    """Lambda^2 element of two letter vectors (dicts) -> dict over pairs"""
    out = {}
    for i, a in u.items():
        for j, b in v.items():
            w = wedge2(i, j)
            if w is None:
                continue
            out[w[1]] = out.get(w[1], 0) + w[0]*a*b
    return {k: x for k, x in out.items() if x != 0}


def outer(X, Y):
    out = {}
    for P, a in X.items():
        for Q, b in Y.items():
            if P == Q:
                continue
            sgn = 1
            PP, QQ = P, Q
            if PP > QQ:
                PP, QQ, sgn = QQ, PP, -1
            out[(PP, QQ)] = out.get((PP, QQ), 0) + sgn*a*b
    return {k: v for k, v in out.items() if v != 0}


def lincomb(*pairs):
    out = {}
    for cf, D in pairs:
        for k, v in D.items():
            out[k] = out.get(k, 0) + cf*v
    return {k: v for k, v in out.items() if v != 0}


if __name__ == '__main__':
    random.seed(7)
    # ---- controls: O kills Li4(f) and shuffle products of random lower-weight symbols
    ok = True
    for _ in range(200):
        f = {i: random.randint(-2, 2) for i in range(3)}
        g = {i: random.randint(-2, 2) for i in range(3)}
        if O(tensor([g, f, f, f])):
            ok = False                                       # symbol of Li4 with 1-f -> g (any letters)
        rnd = lambda L: {tuple(random.randrange(3) for _ in range(L)): Fr(random.randint(-3, 3)) for _ in range(4)}
        if O(shuffle(rnd(1), rnd(3))) or O(shuffle(rnd(2), rnd(2))) or O(shuffle(shuffle(rnd(1), rnd(1)), rnd(2))):
            ok = False
    print('(0) controls: O annihilates Li4-type symbols and all shuffle products (200 random trials each):', ok)
    # check of the symbol rule against two classical identities: Li2(c) = H_{0,1}; Li4(-c) = -H_{0,0,0,-1}
    print('    S(H_{0,1}) =', hsym((0, 1)), ' (expected -[1-c](x)[c] = {(1,0): -1})')
    print('    S(H_{0,0,0,-1}) =', hsym((0, 0, 0, -1)), ' (Li4(-c) = -H_{0,0,0,-1}; S(Li4(-c)) = -[1+c](x)[c]^3)')

    R3 = load('hpl_fit_k3_(0,0,0).json')
    comp = {}
    for e in (-1, 0, 1):
        S = {}
        for ee, w, v in R3:
            if ee == e and len(w) == 4:
                addS(S, hsym(w), v)
        comp[e] = S
    words = {e: sorted({w for ee, w, v in R3 if ee == e and len(w) == 4}) for e in (-1, 0, 1)}
    print('(1) weight-4 words per prefactor:', {e: len(words[e]) for e in words},
          '; distinct overall:', len(set().union(*map(set, words.values()))),
          '; prefactors present at weight 4:', sorted({ee for ee, w, v in R3 if len(w) == 4}))
    names = {(0, 1): '[c]^[1-c]', (0, 2): '[c]^[1+c]', (1, 2): '[1-c]^[1+c]'}
    Ob = {}
    for e in (-1, 0, 1):
        Ob[e] = O(comp[e])
        print(f'    O(A_{e}) =', {f'({names[P]}) ^ ({names[Q]})': str(v) for (P, Q), v in Ob[e].items()} or 0)

    # ---- (2) rho forms (function level, constants dropped)
    def rho(z, omz):
        return vec_wedge(z, omz)
    rho_mt = rho({1: 1, 2: -1}, {2: -1})            # z = -(1-c)/(1+c),  1 - z = 2/(1+c) ~ (1+c)^-1
    rho_c2 = rho({0: 2}, {1: 1, 2: 1})              # z = c^2, 1 - z = (1-c)(1+c)
    rho_c = rho({0: 1}, {1: 1})
    rho_mc = rho({0: 1}, {2: 1})
    t1 = lincomb((Fr(-1, 2), outer(rho_mt, rho_c2)))
    t0 = outer(rho_mt, lincomb((1, rho_c), (-1, rho_mc)))
    print('(2) O(A_-1) == -1/2 rho{-(1-c)/(1+c)} ^ rho{c^2}       :', Ob[-1] == t1)
    print('    O(A_0)  == rho{-(1-c)/(1+c)} ^ (rho{c} - rho{-c})  :', Ob[0] == t0)
    print('    O(A_1)  == 0                                        :', Ob[1] == {})
    keys = sorted(set(Ob[-1]) | set(Ob[0]))
    Mx = sp.Matrix([[Ob[-1].get(k, 0), Ob[0].get(k, 0)] for k in keys])
    print('(5) rank of [O(A_-1), O(A_0)] over Q:', Mx.rank())

    # ---- (3) pullback to tau (function level: constants dropped); tau letters t, 1-t, 1+t, 1+t^2
    TL = {0: 't', 1: '1-t', 2: '1+t', 3: '1+t^2'}
    c_pull = {0: {1: 1, 2: 1, 3: -1}, 1: {0: 2, 3: -1}, 2: {3: -1}}   # c = (1-t)(1+t)/(1+t^2); 1-c = 2t^2/(1+t^2); 1+c = 2/(1+t^2)
    s_pull = {0: {0: 1, 3: -1}, 1: {1: 2, 3: -1}, 2: {2: 2, 3: -1}}   # s = 2t/(1+t^2); 1-s = (1-t)^2/(1+t^2); 1+s = (1+t)^2/(1+t^2)
    t = sp.Symbol('t')
    cexpr = (1 - t**2)/(1 + t**2); sexpr = 2*t/(1 + t**2)
    chk = all(sp.simplify(e1 - e2) == 0 for e1, e2 in [
        (cexpr, (1 - t)*(1 + t)/(1 + t**2)), (1 - cexpr, 2*t**2/(1 + t**2)), (1 + cexpr, 2/(1 + t**2)),
        (1 - sexpr, (1 - t)**2/(1 + t**2)), (1 + sexpr, (1 + t)**2/(1 + t**2))])
    inv = sp.simplify(cexpr.subs(t, (1 - t)/(1 + t)) - sexpr) == 0 and sp.simplify(sexpr.subs(t, (1 - t)/(1 + t)) - cexpr) == 0
    fix = sp.simplify(((1 - t)/(1 + t)).subs(t, sp.sqrt(2) - 1) - (sp.sqrt(2) - 1)) == 0
    print('(3) tau factorizations exact:', chk, '; iota: t -> (1-t)/(1+t) swaps c and s:', inv, '; fixes t = sqrt2-1:', fix)

    def pull(S, mp_):
        out = {}
        for word, cf in S.items():
            addS(out, tensor([mp_[i] for i in word], cf))
        return out
    obs = {}
    for e in (-1, 0, 1):
        for nm, mp_ in (('c', c_pull), ('s', s_pull)):
            obs[(e, nm)] = O(pull(comp[e], mp_))
    print('    nonzero components:', {f'A_{e}({nm})': len(v) for (e, nm), v in obs.items()})
    ks = [(-1, 'c'), (0, 'c'), (-1, 's'), (0, 's')]
    rows = sorted(set().union(*[set(obs[k]) for k in ks]))
    Mx = sp.Matrix([[obs[k].get(r, 0) for k in ks] for r in rows])
    print('    rank of O(A_-1(c)), O(A_0(c)), O(A_-1(s)), O(A_0(s)) in Lambda^2(Lambda^2 <t,1-t,1+t,1+t^2>):', Mx.rank())

    # ---- (4) specialisation at c0 = 1/sqrt2 (constant-retaining letters 2, c, 1-c, 1+c -> <sqrt2, a>)
    spec = {'2': (2, 0), 'c': (-1, 0), '1-c': (-1, -1), '1+c': (-1, 1)}
    r2 = sp.sqrt(2); a_ = 1 + r2; c0 = 1/r2
    okspec = all(sp.simplify(val - r2**p*a_**q) == 0 for val, (p, q) in
                 [(sp.Integer(2), spec['2']), (c0, spec['c']), (1 - c0, spec['1-c']), (1 + c0, spec['1+c'])])
    def vec(d):
        return (sum(e*spec[k][0] for k, e in d.items()), sum(e*spec[k][1] for k, e in d.items()))
    def rho0(z, omz):
        p, q = vec(z), vec(omz)
        return p[0]*q[1] - p[1]*q[0]
    table = {'rho{c}': rho0({'c': 1}, {'1-c': 1}), 'rho{-c}': rho0({'c': 1}, {'1+c': 1}),
             'rho{c^2}': rho0({'c': 2}, {'1-c': 1, '1+c': 1}), 'rho{(1-c)/2}': rho0({'1-c': 1, '2': -1}, {'1+c': 1, '2': -1}),
             'rho{-(1-c)/(1+c)}': rho0({'1-c': 1, '1+c': -1}, {'2': 1, '1+c': -1}),
             'rho{(1-c)/(1+c)}': rho0({'1-c': 1, '1+c': -1}, {'2': 1, 'c': 1, '1+c': -1})}
    print('(4) letters at c0 lie in <sqrt2, 1+sqrt2> (exact):', okspec, '; images as multiples of e = sqrt2 ^ a:', table)
