work = '<project>/research/08_theory/M5a_next_coefficient'
open(f'{work}/w6_symbol.py', 'w').write('''"""Constant-retaining point symbol of W6 and the weight-6 Lie-pairing depth test.

Letters: 0 = X = log 2, 1 = Y = log(1+sqrt2). P = [X,Y].
Classical (Li6 + product) symbols have zero (4,2) and (3,3) Lie parts.
Pairings:
  (3,3): against [X,P] wedge [Y,P]
  (4,2): against A4 wedge P for A4 in {[X,[X,P]], [Y,[X,P]], [Y,[Y,P]]}
Controls: 200 random shuffle products, 100 Li6(S-unit) symbols pair to zero;
Li_{4,2}(b,1), Li_{3,3}(b,1) pair to nonzero; weight-5 code reproduces (84, -38 sqrt2).
"""
import json
import random
import sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent.parent
sys.path.insert(0, str(HERE.parent/'A3_W4_proof'))
sys.path.insert(0, str(HERE.parent/'M4c_w5_closed_form'))
from w4_symbol import *

ZEROQ, ONEQ = Q(0), Q(1)

# ---------------------------------------------------------------- expression builders
POLE = {"0": ZEROQ, "a": A, "-a": -A, "b": B, "-b": -B}

def g_expression_from_omega(omega, i, j):
    factors = (LT,) * i + (MT,) * j
    out = {}
    for name, (p, q) in terms_iter(omega[f"{i},{j}"]):
        z = POLE[name]
        coefficient = Q(F(p), F(q))
        for (w, k), d in log_product(*factors).items():
            key = ((z,) + w, k)
            out[key] = out.get(key, ZEROQ) + coefficient * d
    return clean(out)

def terms_iter(lst):
    return [(nm, pq) for nm, pq in lst]

def w6_expression():
    om = json.loads((HERE/'omega5.json').read_text())
    out = {}
    for ij in [(3, 2), (4, 1), (5, 0)]:
        for key, c in g_expression_from_omega(om, *ij).items():
            out[key] = out.get(key, ZEROQ) + c
    return clean(out)

def symbol(expr):
    out = {}
    for (w, k), c in expr.items():
        s = add(g_symbol(ONE, w, A), scale(g_symbol(ONE, w, B), -1))
        for _ in range(k):
            s = shuffle(s, log_symbol(Q(2)))
        out = add(out, scale(s, c))
    return out

# ---------------------------------------------------------------- Lie words and pairings
def bracket(x, y):
    return add(concat(x, y), scale(concat(y, x), -1))

X, Y = {(0,): ONEQ}, {(1,): ONEQ}
P = bracket(X, Y)
B3 = [bracket(X, P), bracket(Y, P)]                       # [X,P], [Y,P]
B4 = [bracket(X, B3[0]), bracket(Y, B3[0]), bracket(Y, B3[1])]  # [X,[X,P]], [Y,[X,P]], [Y,[Y,P]]
assert all(len(w) == 3 for A_ in B3 for w in A_)
assert all(len(w) == 4 for A_ in B4 for w in A_)

def pair_ab(S, A_lie, B_lie):
    """Contract antisymmetrised deconcat with A wedge B (exact, over Q(sqrt2))."""
    la = len(next(iter(A_lie)))
    lb = len(next(iter(B_lie)))
    tot = ZEROQ
    for w, c in S.items():
        v = A_lie.get(w[:la], ZEROQ) * B_lie.get(w[la:], ZEROQ) \
          - B_lie.get(w[:lb], ZEROQ) * A_lie.get(w[lb:], ZEROQ)
        if v:
            tot += c * v
    return tot

def pair33(S):
    return pair_ab(S, B3[0], B3[1])       # [X,P] wedge [Y,P]

def pair42(S):
    return [pair_ab(S, A4, P) for A4 in B4]

# ---------------------------------------------------------------- controls
def random_word_sym(length, nwords, rng):
    out = {}
    for _ in range(nwords):
        w = tuple(rng.choice((0, 1)) for _ in range(length))
        out[w] = out.get(w, 0) + rng.randint(-5, 5)
    return {w: Q(c) for w, c in out.items() if c}

def shuffle_sym(s1, s2):
    return shuffle(s1, s2)

def li6_symbol_random(rng):
    # Li6(z) symbol = -(1-z) z z z z z with z = +- 2^{j/2} a^m; letters in {X,Y}
    j, m = rng.randint(-4, 4), rng.randint(-3, 3)
    xvec = [F(j, 2), F(m)]       # log z = (j/2) X + m Y
    uvec = [rng.randint(-3, 3), rng.randint(-3, 3)]
    while uvec == [0, 0]:
        uvec = [rng.randint(-3, 3), rng.randint(-3, 3)]
    out = {}
    import itertools
    for bits in itertools.product((0, 1), repeat=6):
        c = uvec[bits[0]]
        for t in range(1, 6):
            c *= xvec[bits[t]]
        if c:
            out[bits] = out.get(bits, ZEROQ) + Q(F(c))
    return clean(out)

def depth2_symbol(m, n, x, y):
    """Symbol of Li_{m,n}(x,y) = G_0(0^{m-1}, 1/x, 0^{n-1}, 1/(x y); 1) (depth 2: sign +)."""
    word = (ZEROQ,)*(m-1) + (ONEQ/x,) + (ZEROQ,)*(n-1) + (ONEQ/(x*y),)
    return g_symbol(ZEROQ, word, ONEQ)

if __name__ == '__main__':
    ban_float()
    self_check()
    report = {}
    # ---- controls: shuffle products and Li6 pair to zero
    rng = random.Random(20261022)
    ok = True
    n_shuf = 0
    for _ in range(200):
        k = rng.choice([1, 2, 3])
        s1 = random_word_sym(k, 3, rng)
        s2 = random_word_sym(6 - k, 3, rng)
        pr = shuffle_sym(s1, s2)
        ok &= (pair33(pr) == ZEROQ) and all(v == ZEROQ for v in pair42(pr))
        n_shuf += 1
    print(f"control: {n_shuf} random shuffle products pair to zero: {ok}", flush=True)
    report['shuffle_products_zero'] = ok
    ok6 = True
    n_li6 = 0
    for _ in range(100):
        s = li6_symbol_random(rng)
        ok6 &= (pair33(s) == ZEROQ) and all(v == ZEROQ for v in pair42(s))
        n_li6 += 1
    print(f"control: {n_li6} Li6(S-unit) symbols pair to zero: {ok6}", flush=True)
    report['li6_symbols_zero'] = ok6
    # ---- depth-2 atoms pair to nonzero
    s42 = depth2_symbol(4, 2, B, ONEQ)
    s33 = depth2_symbol(3, 3, B, ONEQ)
    d42 = [pair33(s42)] + pair42(s42)
    d33 = [pair33(s33)] + pair42(s33)
    print("Li_{4,2}(b,1) pairings (33, 42a,42b,42c):", [str(v) for v in d42], flush=True)
    print("Li_{3,3}(b,1) pairings (33, 42a,42b,42c):", [str(v) for v in d33], flush=True)
    report['li42_b1_pairings'] = [str(v) for v in d42]
    report['li33_b1_pairings'] = [str(v) for v in d33]
    # ---- weight-5 control: reproduce (84, -38 sqrt2)
    import w5_symbol as W5
    s5 = W5.symbol()
    o5 = [pair_ab(s5, W5.L3[0], W5.L2), pair_ab(s5, W5.L3[1], W5.L2)]
    print("weight-5 control pairings ([X,P]xP, [Y,P]xP):", [str(v) for v in o5],
          " (expect 84 and -38 sqrt2)", flush=True)
    report['w5_control'] = [str(v) for v in o5]
    # ---- the W6 symbol and its pairings
    expr = w6_expression()
    print("W6 gpl terms:", len(expr), flush=True)
    S6 = symbol(expr)
    p33 = pair33(S6)
    p42 = pair42(S6)
    print("W6 (3,3) pairing [X,P]x[Y,P]:", str(p33), flush=True)
    print("W6 (4,2) pairings ([X,[X,P]]xP, [Y,[X,P]]xP, [Y,[Y,P]]xP):",
          [str(v) for v in p42], flush=True)
    report['W6_pair33'] = str(p33)
    report['W6_pair42'] = [str(v) for v in p42]
    report['alphabet'] = ['log2', 'log(1+sqrt2)']
    data = dict(report)
    data['symbol'] = {''.join(map(str, w)): [str(c.p), str(c.q)] for w, c in sorted(S6.items())}
    (HERE/'symbol6.json').write_text(json.dumps(data, indent=2) + '\\n')
    print("saved symbol6.json", flush=True)
''')
print("written w6_symbol.py")