"""Independent check of the W6 Lie pairings from symbol6.json (mirrors M4c planner_symbol_pairing.py).

Letters: 0 = X = log2, 1 = Y = log(1+sqrt2). Entries [rational, sqrt2-part].
delta_{a,b} pairing with A (x) B (A in Lie^a, B in Lie^b), a+b=6:
  sum_w s(w) [ <w[:a],A><w[a:],B> - <w[:b],B><w[b:],A> ]  (second term swaps via antisymmetry).
Kills shuffle products and classical Li6. Also redoes all controls and the W5 (84,-38 sqrt2) check.
"""
import json, itertools, random
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent
sym6 = json.loads((HERE / 'symbol6.json').read_text())['symbol']
S6 = {w: (F(p), F(q)) for w, (p, q) in sym6.items()}
sym5 = json.loads((HERE.parent/'M4c_w5_closed_form/symbol.json').read_text())['symbol']
S5 = {w: (F(p), F(q)) for w, (p, q) in sym5.items()}

def lie(expr):
    d = {}
    for c, w in expr: d[w] = d.get(w, 0) + c
    return d
def bracket(A, B):
    out = {}
    for a, ca in A.items():
        for b, cb in B.items():
            out[a+b] = out.get(a+b, 0) + ca*cb
            out[b+a] = out.get(b+a, 0) - ca*cb
    return {w: c for w, c in out.items() if c}
X, Y = lie([(1, '0')]), lie([(1, '1')])
P = bracket(X, Y)
XP, YP = bracket(X, P), bracket(Y, P)
B4 = [bracket(X, XP), bracket(Y, XP), bracket(Y, YP)]

def pair(S, A, B):
    la = len(next(iter(A))); lb = len(next(iter(B)))
    tot = [F(0), F(0)]
    for w, c in S.items():
        v = A.get(w[:la], 0)*B.get(w[la:], 0) - B.get(w[:lb], 0)*A.get(w[lb:], 0)
        if v:
            tot[0] += v*c[0]; tot[1] += v*c[1]
    return tot

def shuffle(a, b):
    if not a: return {b: 1}
    if not b: return {a: 1}
    out = {}
    for w, c in shuffle(a[1:], b).items(): out[a[0]+w] = out.get(a[0]+w, 0) + c
    for w, c in shuffle(a, b[1:]).items(): out[b[0]+w] = out.get(b[0]+w, 0) + c
    return out
def tensor_shuffle(s1, s2):
    out = {}
    for w1, c1 in s1.items():
        for w2, c2 in s2.items():
            for w, c in shuffle(w1, w2).items():
                out[w] = out.get(w, 0) + c*c1*c2
    return out
def to_q(s): return {w: (F(c), F(0)) for w, c in s.items()}

random.seed(7)
ok = True
for _ in range(200):
    k = random.choice([1, 2, 3])
    s1 = {''.join(random.choice('01') for _ in range(k)): random.randint(-5, 5) for _ in range(3)}
    s2 = {''.join(random.choice('01') for _ in range(6-k)): random.randint(-5, 5) for _ in range(3)}
    pr = to_q(tensor_shuffle(s1, s2))
    ok &= pair(pr, XP, YP) == [0, 0] and all(pair(pr, A4, P) == [0, 0] for A4 in B4)
print('independent: 200 shuffle products pair to zero:', ok)
ok = True
for _ in range(100):
    u = (random.randint(-3, 3), random.randint(-3, 3)); x = (random.randint(-3, 3), random.randint(-3, 3))
    s = {}
    for w in itertools.product('01', repeat=6):
        c = u[int(w[0])]*x[int(w[1])]*x[int(w[2])]*x[int(w[3])]*x[int(w[4])]*x[int(w[5])]
        if c: s[''.join(w)] = (F(c), F(0))
    ok &= pair(s, XP, YP) == [0, 0] and all(pair(s, A4, P) == [0, 0] for A4 in B4)
print('independent: 100 Li6 symbols pair to zero:', ok)
print('W6 (3,3) pairing [X,P]x[Y,P]:', pair(S6, XP, YP))
print('W6 (4,2) pairings B4 x P:', [pair(S6, A4, P) for A4 in B4])
print('W5 control pairings ([X,P]xP, [Y,P]xP):', pair(S5, XP, P), pair(S5, YP, P), ' (expect [84,0],[0,-38])')
