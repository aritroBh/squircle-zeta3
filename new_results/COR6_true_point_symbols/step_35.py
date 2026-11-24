from fractions import Fraction as Fr
import importlib, cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *
import random

PRIMES4 = (2, 3, 5, 7)
def vecQ4(x):
    x = Fr(x)
    if x < 0: x = -x
    e = []
    for p in PRIMES4:
        a = 0
        n, d = x.numerator, x.denominator
        while n % p == 0: n //= p; a += 1
        while d % p == 0: d //= p; a -= 1
        e.append(Fr(a))
    n, d = x.numerator, x.denominator
    for p in PRIMES4:
        while n % p == 0: n //= p
        while d % p == 0: d //= p
    assert n == 1 and d == 1, x
    return tuple(e)

def show(S):
    nm = {}
    def nm1(l):
        parts = []
        for p, a in zip(PRIMES4, l):
            if a: parts.append(f"{a}[{p}]")
        return "+".join(parts) if parts else "0"
    return {tuple(nm1(l) for l in k): str(v) for k, v in S.items()}

S = sym_H((0, 0, 0, 1), Fr(1, 2), vecQ4)
print("S(H_{0,0,0,1}(1/2)):", show(S))
vx, vomx = vecQ4(Fr(1, 2)), vecQ4(Fr(1, 2))
print("expect        :", show({(lneg(vomx), vx, vx, vx): Fr(1)}))

# shuffle test, correct multiplicities
ok = True
random.seed(1)
for trial in range(40):
    lu = random.choice([1, 2])
    lv_ = 3 - lu
    u = tuple(random.choice([-1, 0, 1]) for _ in range(lu))
    v = tuple(random.choice([-1, 0, 1]) for _ in range(lv_))
    lhs = sshuffle(sym_H(u, Fr(1, 2), vecQ4), sym_H(v, Fr(1, 2), vecQ4))
    rhs = {}
    for w in shuffles(u, v):
        rhs = sadd(rhs, sym_H(w, Fr(1, 2), vecQ4))
    if lhs != rhs:
        ok = False
        print("FAIL", u, v)
        print(" lhs", show(lhs)); print(" rhs", show(rhs))
        break
print("shuffle consistency:", "PASS" if ok else "FAIL")
