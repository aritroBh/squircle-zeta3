# Planner check of COR5 (independent code):
# (1) the {2,3}-unit equation x + y = 1 by brute force; (2) the symbol of F(1/2) modulo shuffle products equals that of
# -2Li4(-2) + 1/2 Li4(1/3) + 1/2 Li4(2/3) (Lie projector rho, whose kernel is the shuffle ideal); (3) the W3 relation at 300 digits.
import json, itertools
from fractions import Fraction as Fr
import mpmath as mp
# (1)
sols = set()
for a, b, s in itertools.product(range(-12, 13), range(-8, 9), (1, -1)):
    x = s*Fr(2)**a*Fr(3)**b; y = 1 - x
    if y == 0: continue
    n, d = abs(y.numerator), y.denominator
    for p in (2, 3):
        while n % p == 0: n //= p
        while d % p == 0: d //= p
    if n == 1 and d == 1: sols.add(x)
print('(1) #solutions:', len(sols), sorted(sols))
# (2) symbols in V = Q[2] + Q[3]; vectors (e2, e3)
D = json.load(open('<project>/research/10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json'))
c0 = Fr(1, 2)
LET = {0: ((-1, 0), 1), 1: ((-1, 0), -1), -1: ((-1, 1), 1)}   # g(0)=[c], g(1)=-[1-c], g(-1)=[1+c] at c = 1/2
def tens_add(T, word, co):
    # expand multilinear: word = list of vectors
    for idx in itertools.product(range(2), repeat=len(word)):
        coef = co
        for v, i in zip(word, idx): coef *= v[i]
        if coef: T[idx] = T.get(idx, 0) + coef
SF = {}
for key, val in D:
    e, w = key.split(' ', 1); e = int(e[2:]); w = eval(w[1:])
    if len(w) != 4: continue
    word = [LET[a][0] for a in reversed(w)]; sg = 1
    for a in w: sg *= LET[a][1]
    tens_add(SF, word, sg*Fr(val)*c0**e)
def li4(x1mx, x):   # S(Li4(x)) = -[1-x] x [x] x [x] x [x]
    T = {}; tens_add(T, [x1mx, x, x, x], Fr(-1)); return T
SC = {}
for T, co in ((li4((0, 1), (1, 0)), Fr(-2)), (li4((1, -1), (0, -1)), Fr(1, 2)), (li4((0, -1), (1, -1)), Fr(1, 2))):
    for k, v in T.items(): SC[k] = SC.get(k, 0) + co*v
def rho(T):
    # Lie projector: rho(a1..an) = (n-1)/n [rho(a1..a_{n-1}) x a_n - rho(a2..a_n) x a1], rho(a) = a
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
RF, RC = rho(SF), rho(SC)
diff = {k: RF.get(k, 0) - RC.get(k, 0) for k in set(RF) | set(RC)}
print('(2) rho(S(F(1/2))) nonzero terms:', len(RF), '| rho(S(F) - combo) == 0:', all(v == 0 for v in diff.values()))
# control: a shuffle product must project to 0
sh = {}
for w in [((0,), (1,), (0,), (1,))]: pass
ctrl = {}; tens_add(ctrl, [(1, 0), (0, 1)], Fr(1))
def shuffle(u, v):
    if not u: return [v]
    if not v: return [u]
    return [(u[0],) + s for s in shuffle(u[1:], v)] + [(v[0],) + s for s in shuffle(u, v[1:])]
prod = {}
for s in shuffle((0, 1), (1, 1)): prod[s] = prod.get(s, 0) + 1
print('    control rho(shuffle (0,1)sh(1,1)) =', rho(prod), '| rho(Li4(-2)) nonzero:', bool(rho(li4((0, 1), (1, 0)))))
# (3) W3 relation
mp.mp.dps = 320
W3 = mp.mpf(open('<project>/research/08_theory/COR5_corner_classicality/hex_W3.txt').read().split()[0][:330])
l2, l3, z3, pi2 = mp.log(2), mp.log(3), mp.zeta(3), mp.pi**2
rhs = (-mp.mpf(15)/4*z3 - 3*mp.polylog(3, mp.mpf(1)/3) - 8*mp.polylog(3, -mp.mpf(1)/2) - pi2*l2/12 + 4*mp.polylog(2, mp.mpf(1)/3)*l2
       - mp.mpf(3)/2*mp.polylog(2, mp.mpf(1)/3)*l3 + l2**3/3 - mp.mpf(5)/2*l2**2*l3 + 2*l2*l3**2 - l3**3/4)
print('(3) W3 - closed form =', mp.nstr(W3 - rhs, 5))
