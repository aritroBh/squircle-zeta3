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

def normT(S, dim):
    """multilinear normal form"""
    return expand_sym(S, dim)

# variant K: rightmost-outermost, append-right; H_w = (-1)^#1s K(reversed w)
def sym_H_K(word, x, vec):
    s = -1 if sum(1 for a in word if a == 1) % 2 else 1
    return sscale(s, sym_J(tuple(reversed(word)), x, vec, prepend=False, inv=False))

random.seed(1)
print("== shuffle consistency, variant K (reversed, append-right) ==")
ok = True
for trial in range(40):
    lu = random.choice([1, 2])
    lv_ = 3 - lu
    u = tuple(random.choice([-1, 0, 1]) for _ in range(lu))
    v = tuple(random.choice([-1, 0, 1]) for _ in range(lv_))
    lhs = sshuffle(sym_H_K(u, Fr(1, 2), vecQ4), sym_H_K(v, Fr(1, 2), vecQ4))
    rhs = {}
    for w in shuffles(u, v):
        rhs = sadd(rhs, sym_H_K(w, Fr(1, 2), vecQ4))
    if normT(lhs, 4) != normT(rhs, 4):
        ok = False
        print("FAIL", u, v)
print("variant K shuffle:", "PASS" if ok else "FAIL")

print("== shuffle consistency, variant J (prepend-left) ==")
ok = True
for trial in range(40):
    lu = random.choice([1, 2])
    lv_ = 3 - lu
    u = tuple(random.choice([-1, 0, 1]) for _ in range(lu))
    v = tuple(random.choice([-1, 0, 1]) for _ in range(lv_))
    lhs = sshuffle(sym_H(u, Fr(1, 2), vecQ4), sym_H(v, Fr(1, 2), vecQ4))
    rhs = {}
    for w in shuffles(u, v):
        rhs = sadd(rhs, sym_H(w, Fr(1, 2), vecQ4))
    if normT(lhs, 4) != normT(rhs, 4):
        ok = False
        print("FAIL", u, v)
print("variant J shuffle:", "PASS" if ok else "FAIL")

# Li4 tensor check with normal form
print("== Li4 tensor check (multilinear normal form) ==")
for var, fn in [("J", lambda w, x: sym_H(w, x, vecQ4)), ("K", lambda w, x: sym_H_K(w, x))]:
    good = True
    for xs in ["1/2", "1/3", "2/3", "-2", "1/9", "3/4", "5/7", "-3/5"]:
        x = Fr(xs)
        S = fn((0, 0, 0, 1), x)
        vx, vomx = vecQ4(x), vecQ4(1 - x)
        SLi = {(lneg(vomx), vx, vx, vx): Fr(1)}
        if normT(S, 4) != normT(SLi, 4): good = False
    print(f"variant {var}: Li4 tensors all match: {good}")
