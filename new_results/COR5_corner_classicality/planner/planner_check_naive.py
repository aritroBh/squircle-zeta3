# Is the naive specialisation of the function-level HPL symbol at c = 1/2 the true symbol of the number?
# Compare rho(naive S(F(1/2))) with rho(S(Li4 part of COR5's verified closed form)). Products and Li3*log terms have rho = 0.
import itertools
from fractions import Fraction as Fr
exec(open('planner_check.py').read().split('# (3) W3 relation')[0].split('# (2) symbols')[0])  # reuse nothing heavy
import json
D = json.load(open('<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
def tens_add(T, word, co):
    for idx in itertools.product(range(2), repeat=len(word)):
        coef = co
        for v, i in zip(word, idx): coef *= v[i]
        if coef: T[idx] = T.get(idx, 0) + coef
def rho(T):
    memo = {}
    def r(word):
        if word in memo: return memo[word]
        n = len(word)
        if n == 1: out = {word: Fr(1)}
        else:
            out = {}
            for k, v in r(word[:-1]).items(): out[k + (word[-1],)] = out.get(k + (word[-1],), 0) + Fr(n-1, n)*v
            for k, v in r(word[1:]).items(): out[k + (word[0],)] = out.get(k + (word[0],), 0) - Fr(n-1, n)*v
        memo[word] = out; return out
    tot = {}
    for w, co in T.items():
        for k, v in r(w).items(): tot[k] = tot.get(k, 0) + co*v
    return {k: v for k, v in tot.items() if v != 0}
def vec(x):   # x = +-2^a 3^b as Fraction -> (a, b)
    x = abs(x); a = b = 0; n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1; return (a, b)
def li4(x):
    T = {}; tens_add(T, [vec(1 - x), vec(x), vec(x), vec(x)], Fr(-1)); return T
c0 = Fr(1, 2); LET = {0: ((-1, 0), 1), 1: ((-1, 0), -1), -1: ((-1, 1), 1)}
SF = {}
for key, val in D:
    e, w = key.split(' ', 1); e = int(e[2:]); w = eval(w[1:])
    if len(w) != 4: continue
    word = [LET[a][0] for a in reversed(w)]; sg = 1
    for a in w: sg *= LET[a][1]
    tens_add(SF, word, sg*Fr(val)*c0**e)
COMB = {}
for x, co in ((Fr(-2), Fr(-21)), (Fr(1, 3), Fr(-65, 8)), (Fr(2, 3), Fr(13, 2)), (Fr(1, 2), Fr(137, 6)), (Fr(1, 9), Fr(39, 64))):
    for k, v in li4(x).items(): COMB[k] = COMB.get(k, 0) + co*v
rn, rc = rho(SF), rho(COMB)
keys = sorted(set(rn) | set(rc))
print('rho(naive S(F(1/2))) :', {k: str(rn.get(k, 0)) for k in keys})
print('rho(Li4 part of verified closed form):', {k: str(rc.get(k, 0)) for k in keys})
print('equal:', all(rn.get(k, 0) == rc.get(k, 0) for k in keys))
