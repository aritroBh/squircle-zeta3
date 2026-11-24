from fractions import Fraction as Fr
import importlib, cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *
import json

base = "<project>/research/"
D = json.load(open(base + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))

def vecQ(x):
    x = abs(Fr(x)); a = b = 0
    n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1, x
    return (Fr(a), Fr(b))

c0 = Fr(1, 2)
# Omega(1/2): sum over weight-4 terms of coeff * c0^e * S(H_w(1/2))
SOm = {}
for key, val in D:
    e, w = key.split(' ', 1); e = int(e[2:]); w = eval(w[1:])
    if len(w) != 4: continue
    SOm = sadd(SOm, sscale(Fr(val) * c0**e, sym_H(w, c0, vecQ)))

# expand to basis words (0=[2], 1=[3]) and apply planner's rho
SE = expand_sym(SOm, 2)
# convert to index words for rho
def to_idx(SE):
    out = {}
    for w, c in SE.items():
        iw = tuple(0 if l == (Fr(1), Fr(0)) else 1 for l in w)
        out[iw] = out.get(iw, 0) + c
    return out
r = rho(to_idx(SE))
keys = [(1, 0, 0, 0), (1, 1, 0, 0), (1, 1, 1, 0)]
print("rho(S(Omega(1/2))) Lie coords (3222,3322,3332):", tuple(r.get(k, 0) for k in keys))
print("expect (-55/8, 13/8, 13/8)")
print("other keys:", {k: v for k, v in r.items() if k not in keys})

# delta_22 (trivially 0 in dim 2 since Lambda^2 B_2 = 0: B_2 = Lambda^2 V is 1-dim)
d = delta22(SOm, 2)
print("delta_22(Omega(1/2)):", d)

# also whole R3 weight-4 part is Omega; class check vs COR5 closed form Li4 part:
def li4_class(x):
    vx, vomx = vecQ(x), vecQ(1 - x)
    return {(lneg(vomx), vx, vx, vx): Fr(1)}
COMB = {}
for x, co in ((Fr(-2), Fr(-21)), (Fr(1, 3), Fr(-65, 8)), (Fr(2, 3), Fr(13, 2)), (Fr(1, 2), Fr(137, 6)), (Fr(1, 9), Fr(39, 64))):
    COMB = sadd(COMB, sscale(co, li4_class(x)))
rC = rho(to_idx(expand_sym(COMB, 2)))
print("rho(Li4 part of COR5 closed form):", tuple(rC.get(k, 0) for k in keys))
