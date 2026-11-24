from fractions import Fraction as Fr
import importlib, time
import cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe
for m in (cor6_nf, cor6_sunit, cor6_sym, cor6_ft, cor6_pipe): importlib.reload(m)
from cor6_nf import NF, primes_above
from cor6_sunit import exponent_vector
from cor6_sym import *
from cor6_ft import *
import sympy as sp

BASE = "<project>/research/"
import json
D = json.load(open(BASE + "10_frontier/G1_zeta3_any_angle/hpl_fit_k3_(0,0,0).json"))
W4TERMS = [(int(k.split(' ', 1)[0][2:]), eval(k.split(' ', 1)[1][1:]), Fr(v))
           for k, v in D if len(eval(k.split(' ', 1)[1][1:])) == 4]

def letter_span_basis(F, c0, S0):
    """maximal independent subset of {2, c0, 1-c0, 1+c0} spanning all four (verified)."""
    one = F.one()
    cands = [F.from_int(2), c0, F.sub(one, c0), F.add(one, c0)]
    basis = []
    for a in cands:
        if not basis:
            basis.append(a); continue
        if exponent_vector(F, a, basis, dps=120, S0=S0) is not None:
            continue
        basis.append(a)
    return basis

def run_letter_span(label, minpoly, c0, verbose=True):
    F = NF(minpoly)
    one = F.one()
    S0 = set([2])
    for N in (F.norm(c0), F.norm(F.sub(one, c0)), F.norm(F.add(one, c0))):
        S0 |= set(sp.factorint(abs(N.numerator)).keys()) | set(sp.factorint(N.denominator).keys())
    S0 = sorted(p for p in S0 if p > 1)
    basis = letter_span_basis(F, c0, S0)
    letters = {}
    def vec(el):
        if el not in letters:
            ev = exponent_vector(F, el, basis, S0=S0)
            assert ev is not None, ("factor failed", el, label)
            letters[el] = ev
        return letters[el]
    zero = F.zero()
    L = {0: zero, 1: one, -1: F.neg(one)}
    def symH(word):
        wordF = tuple(L[a] for a in word)
        s = -1 if sum(1 for a in word if a == 1) % 2 else 1
        return sscale(s, sym_J(wordF, c0, vec, a0=zero, sub=F.sub, iszero=lambda a: all(x == 0 for x in a)))
    SOm = {}
    for e, w, val in W4TERMS:
        ce = F.pow(c0, e)
        SOm = ft_add(SOm, ft_scale(F, F.scale(val, ce), symH(w)))
    dim = len(basis)
    d22 = ft_delta22(SOm, F.d, dim)
    if verbose:
        print(f"{label}: deg {F.d}, letter-span dim {dim}, delta_22 #components {len(d22)}")
    return F, c0, S0, basis, SOm, d22

# n=7: 2 sin(pi/7): x^6-7x^4+14x^2-7
t0 = time.time()
F7, c7, S07, b7, S7, d22_7 = run_letter_span("c=sin(pi/7) (n=7)", [-7, 0, 14, 0, -7, 0, 1], (Fr(0), Fr(1, 2)) + (Fr(0),)*4)
print(f"time {time.time()-t0:.1f}s, basis: {b7}")
for k, v in sorted(d22_7.items())[:8]:
    print("   ", k, v)
