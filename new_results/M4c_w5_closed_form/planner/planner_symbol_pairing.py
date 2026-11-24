"""Planner check: recompute the (3,2) Lie-pairing obstruction from symbol.json, independently of w5_symbol.py.

Symbol letters: 0 = X = log 2, 1 = Y = log(1+sqrt2). Entries are [rational, sqrt2-part].
delta_{3,2} pairing with A (x) B (A in Lie^3, B in Lie^2):
  sum_w s(w) [ <w[:3],A><w[3:],B> - <w[2:],A><w[:2],B> ]
This kills shuffle products and every classical Li5 (whose symbol ends in x(x)x).
Also checks those two facts on random products and on Li5 symbols.
"""
import json, itertools, random
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
sym = json.loads((HERE / 'symbol.json').read_text())['symbol']
S = {w: (F(p), F(q)) for w, (p, q) in sym.items()}

def lie(expr):  # dict word->coeff from list of (coeff, word)
    d = {}
    for c, w in expr: d[w] = d.get(w, 0) + c
    return d
P = lie([(1, '01'), (-1, '10')])
Q = lie([(1, '001'), (-2, '010'), (1, '100')])          # [X,[X,Y]]
T = lie([(2, '101'), (-1, '110'), (-1, '011')])          # [Y,[X,Y]]

def pair(S, A, B):
    tot = [F(0), F(0)]
    for w, c in S.items():
        v = A.get(w[:3], 0) * B.get(w[3:], 0) - A.get(w[2:], 0) * B.get(w[:2], 0)
        if v:
            tot[0] += v * c[0]; tot[1] += v * c[1]
    return tot

def shuffle(a, b):
    if not a: return {b: 1}
    if not b: return {a: 1}
    out = {}
    for w, c in shuffle(a[1:], b).items(): out[a[0] + w] = out.get(a[0] + w, 0) + c
    for w, c in shuffle(a, b[1:]).items(): out[b[0] + w] = out.get(b[0] + w, 0) + c
    return out

def tensor_shuffle(s1, s2):
    out = {}
    for w1, c1 in s1.items():
        for w2, c2 in s2.items():
            for w, c in shuffle(w1, w2).items():
                out[w] = out.get(w, 0) + c * c1 * c2
    return out

def to_q(s): return {w: (F(c), F(0)) for w, c in s.items()}

random.seed(1)
ok = True
for _ in range(200):  # products of weights (1,4),(2,3)
    k = random.choice([1, 2])
    s1 = {''.join(random.choice('01') for _ in range(k)): random.randint(-5, 5) for _ in range(3)}
    s2 = {''.join(random.choice('01') for _ in range(5 - k)): random.randint(-5, 5) for _ in range(3)}
    pr = to_q(tensor_shuffle(s1, s2))
    ok &= pair(pr, Q, P) == [0, 0] and pair(pr, T, P) == [0, 0]
print('shuffle products pair to zero:', ok)
ok = True
for _ in range(200):  # Li5 symbol (1-x)(x)x(x)x(x)x with letters as vectors over {X,Y}
    u = (random.randint(-3, 3), random.randint(-3, 3)); x = (random.randint(-3, 3), random.randint(-3, 3))
    s = {}
    for w in itertools.product('01', repeat=5):
        c = (u[int(w[0])]) * x[int(w[1])] * x[int(w[2])] * x[int(w[3])] * x[int(w[4])]
        if c: s[''.join(w)] = (F(c), F(0))
    ok &= pair(s, Q, P) == [0, 0] and pair(s, T, P) == [0, 0]
print('classical Li5 symbols pair to zero:', ok)
print('W5 pairing with Q(x)P:', [str(v) for v in pair(S, Q, P)])
print('W5 pairing with T(x)P:', [str(v) for v in pair(S, T, P)])
