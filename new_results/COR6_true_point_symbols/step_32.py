import pickle
from fractions import Fraction as Fr
import importlib, cor6_sym
importlib.reload(cor6_sym)
from cor6_sym import *

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

# (a) controls
print("== (a) single-word controls at x=1/2 ==")
S = sym_H((0, 0, -1, 1), Fr(1, 2), vecQ)
print("H_{0,0,-1,1}(1/2) class:", Q2.coords(S), "expect (-6,0,0) = 6*Li4(-2)")

S = sym_H((0, 0, 1, 1), Fr(1, 2), vecQ)
print("Li_{3,1}(1/2,1) class:", Q2.coords(S), "expect 0")

S = sym_H((0, 1, 0, 1), Fr(1, 2), vecQ)
print("Li_{2,2}(1/2,1) class:", Q2.coords(S), "expect 0")

# H_{0,0,0,1}(x) = Li4(x): symbol must be exactly -[1-x][x]^3
print("== H_{0,0,0,1}(x) = Li4(x) tensor check ==")
for xs in ["1/2", "1/3", "2/3", "-2", "1/9", "3/4", "5/7"]:
    x = Fr(xs)
    S = sym_H((0, 0, 0, 1), x, vecQ)
    vx, vomx = vecQ(x), vecQ(1 - x)
    SLi = {(lneg(vomx), vx, vx, vx): Fr(1)}
    print(f"x={xs}: exact match: {S == SLi}")
