from fractions import Fraction as Fr
import importlib, cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *
import pickle, random

base = "<project>/research/"
wc = pickle.load(open(base + "08_theory/COR5_corner_classicality/word_classes.pkl", "rb"))

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
def sym_H_K(word, x, vec):
    s = -1 if sum(1 for a in word if a == 1) % 2 else 1
    return sscale(s, sym_J(tuple(reversed(word)), x, vec, prepend=False, inv=False))

words24 = sorted(wc.keys())
for name, fn in [("J prepend", lambda w: sym_H(w, Fr(1, 2), vecQ)),
                 ("K reversed-append", lambda w: sym_H_K(w, Fr(1, 2), vecQ))]:
    bad = 0
    for w in words24:
        if Q2.coords(fn(w)) != wc[w][0]:
            bad += 1
    print(f"{name}: class mismatches {bad}/24")

# Do J and K give same tensors? Compare on all 24 words (normal form)
same = 0
for w in words24:
    a = expand_sym(sym_H(w, Fr(1, 2), vecQ), 2)
    b = expand_sym(sym_H_K(w, Fr(1, 2), vecQ), 2)
    if a == b: same += 1
print(f"J and K identical tensors on {same}/24 words")

# example where they differ:
for w in words24:
    a = expand_sym(sym_H(w, Fr(1, 2), vecQ), 2)
    b = expand_sym(sym_H_K(w, Fr(1, 2), vecQ), 2)
    if a != b:
        print("differ on", w)
        break
