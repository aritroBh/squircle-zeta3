"""step_12.py -- FUNCTION-level (3,2) part of the weight-5 part of B4, modulo constants.

The weight-5 part of B4 is the 70 weight-5 terms of R4 (step_09): even-sector combinations
c^e (H_w + sgn H_{-w}), |w| = 5. Its function symbol (constant letters dropped) is
  S5(c) = sum q c^e [ S(H_w) + sgn S(H_{-w}) ]   in Q[c,1/c] (x) T^5(V),
V = span{[c],[1-c],[1+c]} (W1 hpl.py's symbol, with the constant letter '2' projected out).

The (3,2) part is computed in Lie_3(V) (x) Lambda^2(V) (dim 8 x 3 = 24) with the
w5_symbol.py pairing convention: for Q in Lie_3, P in Lambda^2 (as antisymmetric tensors),
  pairing(Q (x) P) = sum_{u in supp Q, v in supp P} c_u d_v (s_{uv} - s_{vu}).
Shuffle products and classical Li5 symbols are killed by this pairing (controls).
Question: is the (3,2) part nonzero AS A FUNCTION OF c?
"""
import sys, json
from fractions import Fraction as Fr
sys.path.insert(0, "<project>/research/10_frontier/W1_corner_family/B_taylor_pi4")
import hpl
import g2_fit as gf
from step_10 import load_terms

# ---- function symbol of the weight-5 part, constant letter (index 0) projected out
def proj_symbol(w):
    """S(H_w) mod constants: dict {word over {1,2,3}: Fr} (basis [c],[1-c],[1+c])."""
    out = {}
    for word, cf in hpl.symbol(tuple(w)).items():
        if any(ix == 0 for ix in word):
            continue
        out[word] = out.get(word, 0) + cf
    return {k: v for k, v in out.items() if v}


def S5_func():
    """returns dict {word: dict {exponent e: Fr}} = Laurent coefficients in c."""
    out = {}
    for kind, w, e, q in load_terms(['hpl_fit_k4_(0,0,0).json']):
        if len(w) != 5:
            continue
        sgn = 1 if (e + gf.nz(w)) % 2 == 0 else -1
        for sw, mul in ((w, 1), (gf.neg(w), sgn)):
            for word, cf in proj_symbol(sw).items():
                out.setdefault(word, {})
                out[word][e] = out[word].get(e, 0) + q*cf*mul
    return {k: {e: v for e, v in d.items() if v} for k, d in out.items() if any(d.values())}


# ---- Lie machinery on tensor words over {1,2,3}
def concat(x, y):
    return {x + y: Fr(1)}

def bracket(x, y):
    out = {x + y: Fr(1)}
    out[y + x] = out.get(y + x, 0) - 1
    return {k: v for k, v in out.items() if v}

def br(A, B):
    """[A, B] = AB - BA for tensor elements A, B (dicts)."""
    out = {}
    for u, cu in A.items():
        for v, cv in B.items():
            out[u + v] = out.get(u + v, 0) + cu*cv
            out[v + u] = out.get(v + u, 0) - cu*cv
    return {k: v for k, v in out.items() if v}


X = {(1,): Fr(1)}
Y = {(2,): Fr(1)}
Z = {(3,): Fr(1)}
LIE3 = [br(X, br(X, Y)), br(Y, br(X, Y)),
        br(X, br(X, Z)), br(Z, br(X, Z)),
        br(Y, br(Y, Z)), br(Z, br(Y, Z)),
        br(X, br(Y, Z)), br(Y, br(X, Z))]
WEDGE2 = [br(X, Y), br(X, Z), br(Y, Z)]
LIE3_NAMES = ['[X,[X,Y]]', '[Y,[X,Y]]', '[X,[X,Z]]', '[Z,[X,Z]]',
              '[Y,[Y,Z]]', '[Z,[Y,Z]]', '[X,[Y,Z]]', '[Y,[X,Z]]']
WEDGE2_NAMES = ['X^Y', 'X^Z', 'Y^Z']


def pairing(S5, Q, P):
    """sum c_u d_v (s_{uv} - s_{vu}); result: Laurent polynomial dict {e: Fr}."""
    out = {}
    for u, cu in Q.items():
        for v, dv in P.items():
            for word, sgn in ((u + v, 1), (v + u, -1)):
                for e, cf in S5.get(word, {}).items():
                    out[e] = out.get(e, 0) + sgn*cu*dv*cf
    return {e: v for e, v in out.items() if v}


def fmt_laurent(d):
    if not d:
        return '0'
    return ' + '.join(f'({v})c^{e}' for e, v in sorted(d.items()))


if __name__ == '__main__':
    S5 = S5_func()
    print(f'weight-5 function symbol: {len(S5)} tensor words')
    # controls: shuffle products and Li5 symbols must be killed (spot check at the function level)
    def shuffle_words(u, v):
        if not u: return [v]
        if not v: return [u]
        return [(u[0],) + w for w in shuffle_words(u[1:], v)] + [(v[0],) + w for w in shuffle_words(u, v[1:])]
    import random
    random.seed(3)
    ctl = 0
    for t in range(200):
        lu = random.choice([1, 2, 3, 4])
        lv = 5 - lu
        u = tuple(random.choice([1, 2, 3]) for _ in range(lu))
        v = tuple(random.choice([1, 2, 3]) for _ in range(lv))
        Ssh = {}
        for w in shuffle_words(u, v):
            Ssh.setdefault(w, {0: Fr(0)})
            Ssh[w][0] += 1
        for Q in random.sample(LIE3, 2):
            for P in random.sample(WEDGE2, 1):
                assert not pairing(Ssh, Q, P), 'shuffle control failed'
        ctl += 1
    print(f'shuffle controls passed: {ctl} random products x 2 pairings')
    print()
    print('(3,2) pairing matrix (Laurent polynomials in c):')
    any_nonzero = False
    M = {}
    for i, Q in enumerate(LIE3):
        for j, P in enumerate(WEDGE2):
            r = pairing(S5, Q, P)
            M[(i, j)] = r
            if r:
                any_nonzero = True
                print(f'  {LIE3_NAMES[i]:>10} (x) {WEDGE2_NAMES[j]:>3} : {fmt_laurent(r)}')
    print()
    print('FUNCTION-level (3,2) part of B4 is NONZERO as a function of c:' , any_nonzero)
    json.dump({f'{LIE3_NAMES[i]}|{WEDGE2_NAMES[j]}': {str(e): str(v) for e, v in r.items()}
               for (i, j), r in M.items()}, open('func_delta32.json', 'w'))
