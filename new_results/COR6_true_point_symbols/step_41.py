from fractions import Fraction as Fr
import importlib, cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *
import json, itertools

base = "<project>/research/"
D = json.load(open(base + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))

# planner's naive specialization
def tens_add(T, word, co):
    for idx in itertools.product(range(2), repeat=len(word)):
        coef = co
        for v, i in zip(word, idx): coef *= v[i]
        if coef: T[idx] = T.get(idx, 0) + coef
def vec(x):
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
rn = rho(SF)
print("naive rho full dict:", {k: str(v) for k, v in sorted(rn.items())})
print("reading (0,1,0,0),(0,1,0,1),(0,1,1,0):", (rn.get((0,1,0,0),0), rn.get((0,1,0,1),0), rn.get((0,1,1,0),0)))
print("task says naive = (-5/8, 1/8, 1/8)")

# true symbol reading:
def vecQ(x):
    x = abs(Fr(x)); a = b = 0
    n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1, x
    return (Fr(a), Fr(b))
SOm = {}
for key, val in D:
    e, w = key.split(' ', 1); e = int(e[2:]); w = eval(w[1:])
    if len(w) != 4: continue
    SOm = sadd(SOm, sscale(Fr(val) * c0**e, sym_H(w, c0, vecQ)))
SE = expand_sym(SOm, 2)
T = {}
for w, c in SE.items():
    iw = tuple(0 if l == (Fr(1), Fr(0)) else 1 for l in w)
    T[iw] = T.get(iw, 0) + c
rt = rho(T)
print("true rho full dict:", {k: str(v) for k, v in sorted(rt.items())})
print("reading:", (rt.get((0,1,0,0),0), rt.get((0,1,0,1),0), rt.get((0,1,1,0),0)))
print("task says true = (-55/8, 13/8, 13/8)")
