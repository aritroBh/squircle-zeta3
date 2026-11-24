"""b2b2_test.py -- EXACT (symbol level).
Classicality test for the weight-4 part W_4(c) of the corner coefficient M_3(alpha) as a FUNCTION of c.

W_4(c) = c^{-1} A_{-1}(c) + A_0(c) + c A_1(c)  with A_e pure weight-4 combinations of HPLs H_w(c)
(letters 0,+-1) plus z2*(weight 2) + z3*(weight 1) (zero symbol).  By linear independence of HPLs over Q(c)
(Deneufchatel-Duchamp-Minh-Solomon), a classical expression F(c) = sum p_i(c) Li_4(x_i(c)) + products
with p_i in Q(c) can equal W_4(c) only if every A_e is, modulo products, a Q-combination of Li_4's.
A necessary condition (Goncharov; GSVV 2010) is that the Lambda^2 B_2 component of the cobracket vanishes:
   delta_{2,2}(S) = sum coef * (a^b) ^ (c^d)   for S = sum coef a(x)b(x)c(x)d,
an element of Lambda^2(Lambda^2 A), A = <2, c, 1-c, 1+c>_Q (dimension 15).
Products and Li_4's map to 0 (checked below on test cases).
"""
import itertools, json, re, sys
from fractions import Fraction as Fr
sys.path.insert(0, '.')
import hpl

G1 = 'research/10_frontier/G1_zeta3_any_angle/'
LET = hpl.BASIS
PAIRS = [(i, j) for i in range(4) for j in range(i + 1, 4)]
PIDX = {p: k for k, p in enumerate(PAIRS)}


def wedge(i, j):
    """a_i ^ a_j as (sign, pair index) or None"""
    if i == j:
        return None
    if i < j:
        return (1, PIDX[(i, j)])
    return (-1, PIDX[(j, i)])


def delta22(S):
    out = {}
    for (a, b, c, d), cf in S.items():
        w1 = wedge(a, b); w2 = wedge(c, d)
        if w1 is None or w2 is None:
            continue
        s = w1[0]*w2[0]
        P, Q = w1[1], w2[1]
        if P == Q:
            continue
        if P > Q:
            P, Q, s = Q, P, -s
        out[(P, Q)] = out.get((P, Q), 0) + s*cf
        if out[(P, Q)] == 0:
            del out[(P, Q)]
    return out


def shuffle(S1, S2):
    out = {}
    for w1, c1 in S1.items():
        for w2, c2 in S2.items():
            n1, n2 = len(w1), len(w2)
            for pos in itertools.combinations(range(n1 + n2), n1):
                word = [None]*(n1 + n2)
                it1 = iter(w1); it2 = iter(w2)
                ps = set(pos)
                for k in range(n1 + n2):
                    word[k] = next(it1) if k in ps else next(it2)
                t = tuple(word)
                out[t] = out.get(t, 0) + c1*c2
    return {k: v for k, v in out.items() if v != 0}


def lin(*pairs):
    out = {}
    for cf, S in pairs:
        for k, v in S.items():
            out[k] = out.get(k, 0) + cf*v
    return {k: v for k, v in out.items() if v != 0}


def load(name):
    out = []
    for nm, v in json.load(open(G1 + name)):
        m = re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
        if not m:
            continue
        out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), Fr(v)))
    return out


def fmt(D):
    names = [f'({LET[i]})^({LET[j]})' for (i, j) in PAIRS]
    return '  ' + '\n  '.join(f'{v} * [{names[P]}] ^ [{names[Q]}]' for (P, Q), v in sorted(D.items()))


if __name__ == '__main__':
    # --- sanity: products and classical Li4 give 0
    S_li2x = {(2, 1): Fr(-1)}           # Li2(x)  = -[1-x](x)[x]
    S_li2mx = {(3, 1): Fr(-1)}          # Li2(-x) = -[1+x](x)[x]
    print('delta22( Li2(x)*Li2(-x) ) =', delta22(shuffle(S_li2x, S_li2mx)))
    S_h1m1 = hpl.symbol((1, -1))
    print('delta22( H_{1,-1}*H_{0,-1} ) =', delta22(shuffle(S_h1m1, hpl.symbol((0, -1)))))
    print('delta22( log(1+x)*H_{0,1,-1} ) =', delta22(shuffle({(3,): Fr(1)}, hpl.symbol((0, 1, -1)))))
    # Li4((1-x)/2): -[1-z](x)[z]^3 with z = (1-x)/2, 1-z = (1+x)/2 : letters as vectors
    def vec_word(vs, cf):
        out = {}
        for idx in itertools.product(*[[(k, e) for k, e in enumerate(v) if e] for v in vs]):
            key = tuple(k for k, _ in idx)
            val = cf
            for _, e in idx:
                val *= e
            out[key] = out.get(key, 0) + val
        return {k: v for k, v in out.items() if v != 0}
    z = (-1, 0, 1, 0); omz = (-1, 0, 0, 1)
    print('delta22( Li4((1-x)/2) ) =', delta22(vec_word([omz, z, z, z], Fr(-1))))

    # --- the actual test
    R3 = load('hpl_fit_k3_(0,0,0).json')
    total = {}
    for e in (-1, 0, 1):
        S = {}
        for ee, w, v in R3:
            if ee == e and len(w) == 4:
                S = lin((1, S), (v, hpl.symbol(w)))
        D = delta22(S)
        print(f'\nprefactor c^{e}: weight-4 symbol has {len(S)} terms; delta_22 has {len(D)} nonzero components')
        if D:
            print(fmt(D))
    # also the even/odd structure: show W4 total symbol size
