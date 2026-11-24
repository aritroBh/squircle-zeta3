import pickle
from fractions import Fraction as Fr
import importlib, cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *

base = "<project>/research/"
wc = pickle.load(open(base + "08_theory/COR5_corner_classicality/word_classes.pkl", "rb"))

def vecQ(x):
    # x: nonzero Fraction; vector (a,b) with |x| = 2^a 3^b
    x = abs(Fr(x)); a = b = 0
    n, d = x.numerator, x.denominator
    while n % 2 == 0: n //= 2; a += 1
    while d % 2 == 0: d //= 2; a -= 1
    while n % 3 == 0: n //= 3; b += 1
    while d % 3 == 0: d //= 3; b -= 1
    assert n == 1 and d == 1, x
    return (Fr(a), Fr(b))

Q2 = W4Quotient(2)
print("dim Q4:", len(Q2.free), "free words:", Q2.free_words)
c0 = Fr(1, 2)
words24 = sorted(wc.keys())
for prepend, inv in [(True, False), (True, True), (False, False), (False, True)]:
    bad = []
    for w in words24:
        S = sym_H(w, c0, vecQ, prepend=prepend, inv=inv)
        cl = Q2.coords(S)
        if cl != wc[w][0]:  # numeric (true) class
            bad.append((w, cl, wc[w][0]))
    print(f"prepend={prepend} inv={inv}: mismatches {len(bad)}/24")
    for b in bad[:5]: print("   ", b)
