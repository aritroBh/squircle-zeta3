from fractions import Fraction as Fr
import importlib, json
import cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *

BASE = "<project>/research/"
D = json.load(open(BASE + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
W4TERMS = [(int(k.split(' ', 1)[0][2:]), eval(k.split(' ', 1)[1][1:]), Fr(v))
           for k, v in D if len(eval(k.split(' ', 1)[1][1:])) == 4]

def vecQ(x):
    x = abs(Fr(x)); a = b = 0
    n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1, x
    return (Fr(a), Fr(b))

Q2 = W4Quotient(2)
def li4_sym(x):
    vx, vomx = vecQ(x), vecQ(1 - x)
    return {(lneg(vomx), vx, vx, vx): Fr(1)}

SOm3 = {}
for e, w, val in W4TERMS:
    SOm3 = sadd(SOm3, sscale(val * Fr(1, 3)**e, sym_H(w, Fr(1, 3), vecQ)))
clsOm3 = Q2.coords(SOm3)

comb = [(Fr(-8), Fr(-11, 81)), (Fr(-3), Fr(-2)), (Fr(2, 3), Fr(20, 3))]
Scomb = {}
for x, co in comb:
    Scomb = sadd(Scomb, sscale(co, li4_sym(x)))
clsComb = Q2.coords(Scomb)
print("class(Omega(1/3)):", clsOm3)
print("class(combination):", clsComb)
print("MATCH:", clsOm3 == clsComb)
