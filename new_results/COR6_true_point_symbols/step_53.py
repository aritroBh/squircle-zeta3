from fractions import Fraction as Fr
import importlib, json
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft): importlib.reload(m)
from cor6_nf import NF
from cor6_sunit import build_basis, exponent_vector
from cor6_sym import *
from cor6_ft import *

base = "<project>/research/"
D = json.load(open(base + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))

F = NF([-2, 0, 1])
d = 2
basis, pool = build_basis(F, [2], 2, verbose=False)
sq2 = F.theta()
c0 = (Fr(0), Fr(1, 2))
zero, one = F.zero(), F.one()
L0, L1, Lm1 = zero, one, F.neg(one)

letters_needed = {}
def vecF(el):
    if el not in letters_needed:
        ev = exponent_vector(F, el, basis)
        assert ev is not None, ("factor failed", el)
        letters_needed[el] = ev
    return letters_needed[el]
sub = lambda a, b: F.sub(a, b)
iszero = lambda a: all(x == 0 for x in a)
def sym_H_F(word, c0):
    wordF = tuple(L1 if a == 1 else (Lm1 if a == -1 else L0) for a in word)
    s = -1 if sum(1 for a in word if a == 1) % 2 else 1
    return sscale(s, sym_J(wordF, c0, vecF, a0=zero, sub=sub, iszero=iszero))

SOm = {}
for key, val in D:
    e, w = key.split(' ', 1); e = int(e[2:]); w = eval(w[1:])
    if len(w) != 4: continue
    ce = F.pow(c0, e)
    SOm = ft_add(SOm, ft_scale(F, F.scale(Fr(val), ce), sym_H_F(w, c0)))

def li4_ft(q, x):
    vx, vomx = vecF(x), vecF(F.sub(one, x))
    return {(lneg(vomx), vx, vx, vx): q}

b = F.sub(sq2, one)  # sqrt2 - 1
mb = F.neg(b)
mb2 = F.neg(F.mul(b, b))
x1 = c0  # 1/sqrt2
half = F.scale(Fr(1, 2), one)
x3 = F.sub(one, c0)  # 1 - 1/sqrt2
x4 = F.sub(F.from_int(2), sq2)  # 2 - sqrt2
def cn(p, q): return (Fr(p), Fr(q))
RLi4 = ft_add(
    li4_ft(cn(1, -573), b),
    li4_ft(cn(-1, -573), mb),
    li4_ft(cn(0, 56), mb2),
    li4_ft(cn(0, 17), x1),
    li4_ft(cn(Fr(-37, 4), Fr(-17, 16)), half),
    li4_ft(cn(10, 86), x3),
    li4_ft(cn(-10, 86), x4),
)
Q2 = W4Quotient(2)
cOm = ft_coords(Q2, SOm, d)
cR = ft_coords(Q2, RLi4, d)
print("class(Omega(1/sqrt2)):", cOm)
print("class(Li4 part of R) :", cR)
print("MATCH:", cOm == cR)
