"""g8_mutants.py -- the all-n identity is not vacuous: each single change of a coefficient or of an exponent in the data of
Theorem 20 / Theorem 21 makes the comparison in g8_allN.py fail.  (A test of the test; not part of the proof.)"""
import sys
from fractions import Fraction as Fr
import g8_allN as A
from g8_zsum import K, nn, lin, add_into

def compare(k, zm, terms, extra):
    P = A.mu_side(k)
    ev = A.hpl_side(terms, 0)
    for w, c in extra.items(): add_into(ev, w, c)
    diff = lin((1, ev), (-1, P[zm]))
    odd = A.hpl_side(terms, 1)
    lo = A.low_order(terms, lambda N: Fr(0), lo=-4, hi=1)
    return (len(diff) == 0 and len(odd) == 0 and not lo), len(diff), len(odd), len(lo)

base = A.m2_terms_from_statement()
extra = {(): -K(1)/(2*nn**2)}
print("k=2 unmodified:", compare(2, (0,), base, extra))
# mutant a: Y: coefficient 1/2 of c^-1 H_{1,0,1} -> 1
m = [(kd, w, e, (q*2 if (w, e) == ((1, 0, 1), -1) else q)) for kd, w, e, q in base]
print("k=2 mutant a (coefficient of c^-1 H_101 doubled):", compare(2, (0,), m, extra))
# mutant b: a letter changed
m = [(kd, ((0, 1, 1) if (w, e) == ((0, 1, -1), 0) else w), e, q) for kd, w, e, q in base]
print("k=2 mutant b (H_{0,1,-1} replaced by H_{0,1,1}):", compare(2, (0,), m, extra))
# mutant c: drop -1/2 Li_2(c^2)
print("k=2 mutant c (Li_2(c^2) term dropped):", compare(2, (0,), base, {}))
# mutant d: prefactor exponent
m = [(kd, w, (0 if (w, e) == ((1, 0, -1), -1) else e), q) for kd, w, e, q in base]
print("k=2 mutant d (c^-1 H_{1,0,-1} -> c^0 H_{1,0,-1}):", compare(2, (0,), m, extra))
t4 = A.parse_terms('g8_terms_k4_part1.json')
print("k=4 part 1 unmodified:", compare(4, (0, 0, 0), t4, {}))
m = [(kd, w, e, (q + 1 if i == 3 else q)) for i, (kd, w, e, q) in enumerate(t4)]
print("k=4 part 1 mutant (coefficient of term 4 increased by 1):", compare(4, (0, 0, 0), m, {}))
t4z = A.parse_terms('g8_terms_k4_partz2.json')
m = [(kd, w, e, (q*2 if i == 5 else q)) for i, (kd, w, e, q) in enumerate(t4z)]
print("k=4 part zeta(2) mutant (coefficient of term 6 doubled):", compare(4, (1, 0, 0), m, {}))
