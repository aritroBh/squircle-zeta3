import pickle
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

c0 = Fr(1, 2)
for w in [(0,), (1,), (0, 0, 0, 1), (0, 0, 0, -1), (0, 0, -1, 1)]:
    S = sym_J(w, c0, vecQ)
    print("J", w, {tuple(tuple(str(y) for y in l) for l in k): str(v) for k, v in S.items()})
