
import pickle, sys, json
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
import itertools
c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
Z3 = (Fr(0),)*3

def factor_vec(rat):
    num, den = sp.fraction(sp.together(rat))
    vec = [0, 0, 0, 0]
    for poly, sgn in [(sp.expand(num), +1), (sp.expand(den), -1)]:
        if poly == 0: return None
        fac = sp.factor_list(poly)
        base = sp.Rational(fac[0]); b = abs(base)
        if b != 1:
            bn, bd = b.numerator, b.denominator
            if bn & (bn - 1) or bd & (bd - 1): return None
            vec[0] += sgn * (bn.bit_length() - bd.bit_length())
        for f, e in fac[1]:
            f = sp.expand(f)
            if f == c: vec[1] += sgn * e
            elif f == c - 1: vec[2] -= sgn * e
            elif f == c + 1: vec[3] += sgn * e
            elif f.is_Number:
                ff = abs(sp.Rational(f)); fn, fd = ff.numerator, ff.denominator
                if fn & (fn - 1) or fd & (fd - 1): return None
                vec[0] += sgn * e * (fn.bit_length() - fd.bit_length())
            else: return None
    return tuple(Fr(v) for v in vec[1:])

def logfacs(arg):
    """[arg] and [1-arg] as letter vectors"""
    return factor_vec(arg), factor_vec(1 - arg)

def term_symbol(kind, data):
    if kind == 'Li4':
        rv, omrv = logfacs(data)
        return li_sym(4, rv, omrv)
    S = None
    for k, arg in data:
        if k == 'log':
            v = factor_vec(arg)
            s = log_sym(v)
        else:
            rv, omrv = logfacs(arg)
            s = li_sym(int(k[2]), rv, omrv)
        S = s if S is None else sshuffle(S, s)
    return S

bwords = list(itertools.product(range(3), repeat=4))
bidx = {w:i for i,w in enumerate(bwords)}
def expand_letter(l):
    return [(j, l[j]) for j in range(3) if l[j]]
def to81(S, e=0):
    v = [Fr(0)]*81
    for w, cf in S.items():
        coeff = cf.get(e)
        if coeff is None: continue
        terms = [((), Fr(1))]
        for l in w:
            opts = expand_letter(l)
            terms = [(p+(j,), c0*m) for p, c0 in terms for j, m in opts]
        for pp, c0 in terms:
            v[bidx[pp]] += coeff*c0
    return v

terms = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/identity_terms.pkl','rb'))
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81.pkl','rb'))
Rv = D['Rv']
# reconstruct: sum over parsed terms of cf[E]*to81(symbol) == Rv[E] ?
acc = {E: [Fr(0)]*81 for E in (-1,0,1)}
for kind, data, cf in terms:
    S = term_symbol(kind, data)
    for E, v0 in cf.items():
        v = to81(S)
        acc[E] = [x + v0*y for x, y in zip(acc[E], v)]
for E in (-1,0,1):
    diff = max(abs(x - y) for x, y in zip(Rv[E], acc[E]))
    print(f"E={E}: parsed-terms symbol residual:", diff)
