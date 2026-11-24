"""Planner: independent depth test for W6 (weight-6 part of M5), from M5a's gpl_w6.tsv.

Symbol: A3_W4_proof/w4_symbol.py engine (constant-retaining point symbol; letters reduced to X = log 2 (index 0),
Y = log(1+sqrt2) (index 1)). Pairings (my own code):
  (3,3): <S, A (x) B - B (x) A> with A = [X,P], B = [Y,P], P = [X,Y]          (Lambda^2 B3, 1-dim)
  (4,2): <S, C (x) P - P (x) C> for C in {[X,[X,P]], [X,[Y,P]], [Y,[Y,P]]}  (B4 (x) B2, 3-dim)
Classical (Li6 + products) symbols give 0 for all four. Controls: shuffles, Li6, Li_{4,2}, Li_{3,3}, and the weight-5
version of the same code on gpl_w5.tsv must give M4c's (3,2) obstruction (84, -38 sqrt2).
"""
import random, re, sys
from fractions import Fraction as F
from pathlib import Path
TH = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(TH / 'A3_W4_proof'))
from w4_symbol import Q, R, A, B, ZERO, ONE, add, scale, shuffle, log_symbol, concat, g_symbol, li_symbol

def fld(s):
    s = s.strip().strip('()').replace('sqrt2', 's'); v = ZERO
    for t in re.findall(r'[+-]?[^+-]+', s):
        if 's' in t:
            c = t.replace('*s', '').replace('s', ''); c = c if c not in ('', '+', '-') else c + '1'
            v = v + Q(F(c)) * R
        else:
            v = v + Q(F(t))
    return v

def load(path):
    terms = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('#') or line.startswith('coefficient'): continue
        p, q, k, word = line.split('\t')
        w = tuple(fld(x) for x in re.findall(r'\((?:[^()]|\([^()]*\))*\)|[^,() ]+', word.strip()[1:-1]))
        terms.append((Q(F(p)) + Q(F(q)) * R, int(k), w))
    return terms

def symbol_of(terms):
    out = {}
    for c, k, w in terms:
        s = add(g_symbol(ONE, w, A), scale(g_symbol(ONE, w, B), -1))
        for _ in range(k): s = shuffle(s, log_symbol(Q(2)))
        out = add(out, scale(s, c))
    return out

X, Y = {(0,): ONE}, {(1,): ONE}
def br(a, b): return add(concat(a, b), scale(concat(b, a), -1))
P = br(X, Y)
L3 = [br(X, P), br(Y, P)]
L4 = [br(X, br(X, P)), br(X, br(Y, P)), br(Y, br(Y, P))]

def pair(S, U, V):  # <S, U (x) V - V (x) U>
    tot = ZERO
    for u, a in U.items():
        for v, b in V.items():
            tot = tot + a * b * (S.get(u + v, ZERO) - S.get(v + u, ZERO))
    return tot

def tests6(S): return [pair(S, L3[0], L3[1])] + [pair(S, C, P) for C in L4]
def fmt(v): return '(%s, %s)' % (v.p, v.q) if hasattr(v, 'p') else str(v)

# weight-5 control on gpl_w5.tsv: (3,2) pairing against [X,P](x)P and [Y,P](x)P
S5 = symbol_of(load(TH / 'M4a_decomposition' / 'gpl_w5.tsv'))
print('W5 control (3,2):', [fmt(pair(S5, A3, P)) for A3 in L3], ' expect (84,0) and (0,-38)')

random.seed(1)
def rnd(n): return {tuple(random.randint(0, 1) for _ in range(n)): Q(F(random.randint(-5, 5))) for _ in range(4)}
bad = 0
for _ in range(200):
    a = random.randint(1, 5); s = shuffle(rnd(a), rnd(6 - a))
    bad += any(v != ZERO for v in tests6(s))
print('shuffle controls nonzero:', bad, '/ 200')
units = [R - 1, 1 - R, 3 - 2 * R, -3 + 2 * R, R / 2, Q(F(1, 2)), 1 - R / 2, 2 - R, Q(F(1, 2)) - R / 4, -2 + 2 * R]
print('Li6 controls nonzero:', sum(any(v != ZERO for v in tests6(li_symbol(6, z))) for z in units), '/', len(units))
a = A  # 1/b
print('Li_{4,2}(b,1) pairing:', [fmt(v) for v in tests6(g_symbol(ZERO, (ZERO, ZERO, ZERO, a, ZERO, a), ONE))])
print('Li_{3,3}(b,1) pairing:', [fmt(v) for v in tests6(g_symbol(ZERO, (ZERO, ZERO, a, ZERO, ZERO, a), ONE))])

S6 = symbol_of(load(TH / 'M5a_next_coefficient' / 'gpl_w6.tsv'))
print('W6 symbol: %d words' % len(S6))
r = tests6(S6)
print('W6 (3,3) [X,P]^[Y,P]:', fmt(r[0]))
for name, v in zip(['[X,[X,P]]', '[X,[Y,P]]', '[Y,[Y,P]]'], r[1:]): print('W6 (4,2)', name, '(x) P:', fmt(v))
