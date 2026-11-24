from fractions import Fraction as Fr
import importlib, cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *

def vecQ(x):
    x = abs(Fr(x)); a = b = 0
    n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1, x
    return (Fr(a), Fr(b))

def show(S):
    nm = {(1, 0): "[2]", (0, 1): "[3]", (-1, 0): "-[2]", (0, -1): "-[3]",
          (1, -1): "[2]-[3]", (-1, 1): "[3]-[2]", (-2, 0): "-2[2]", (0, 2): "2[3]",
          (2, 0): "2[2]", (0, -2): "-2[3]"}
    return {tuple(nm.get(tuple(l), str(l)) for l in k): str(v) for k, v in S.items()}

S = sym_H((0, 0, 0, 1), Fr(1, 2), vecQ)
print("S(H_{0,0,0,1}(1/2)) =", show(S))
vx, vomx = vecQ(Fr(1, 2)), vecQ(Fr(1, 2))
print("expected:", show({(lneg(vomx), vx, vx, vx): Fr(1)}))
