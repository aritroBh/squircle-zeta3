
import json, pickle, sys
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
            elif f == c - 1: vec[2] += sgn * e
            elif f == c + 1: vec[3] += sgn * e
            elif f.is_Number:
                ff = abs(sp.Rational(f)); fn, fd = ff.numerator, ff.denominator
                if fn & (fn - 1) or fd & (fd - 1): return None
                vec[0] += sgn * e * (fn.bit_length() - fd.bit_length())
            else: return None
    return tuple(Fr(v) for v in vec)

# re-enumerate the 40 S-unit args (same ranges as step13)
args = {}
for a in range(-2, 3):
    for b in range(-3, 4):
        for d in range(-3, 4):
            for e in range(-3, 4):
                for sgn in [1, -1]:
                    r = sgn * sp.Integer(2)**a * c**b * (1 - c)**d * (1 + c)**e
                    v_r = (Fr(a), Fr(b), Fr(d), Fr(e))
                    v_om = factor_vec(sp.expand(1 - r))
                    if v_om is not None:
                        args[(v_r, v_om)] = r
print(len(args), "S-unit args (correct enumeration)")
sa = [(sp.sstr(r), v_r[1:], v_om[1:]) for (v_r, v_om), r in sorted(args.items(), key=lambda kv: str(kv[1]))]
nonconst = [(s, rv, omv) for s, rv, omv in sa if rv != Z3]
print(len(nonconst), "nonconstant")

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

def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)
# S(Omega)
FS = json.load(open(C4 + "F_symbol.json"))
def parse_letter(s):
    digs = []
    i = 0
    while i < len(s):
        if s[i] == '-':
            digs.append(int(s[i:i+2])); i += 2
        else:
            digs.append(int(s[i])); i += 1
    assert digs[0] == 0
    return (Fr(digs[1]), Fr(digs[2]), Fr(digs[3]))
OM = {}
for key, val in FS.items():
    e = int(key.split(" ")[0][2:])
    word = tuple(parse_letter(b) for b in key.split(" ")[1].split("|"))
    OM[word] = padd(OM.get(word, {}), {e: Fr(val)})
OMv = {E: to81(OM, E) for E in (-1,0,1)}
# DT = (c^-1-1)Li22(1/2,1-c) - (c^-1+1)Li22(1/2,1+c)
def drop_zero(S):
    return {w: cf for w, cf in S.items() if all(any(x for x in l) for l in w)}
S22m = drop_zero(li22_sym(Z3, mv, Z3, cv, pv))
S22p = drop_zero(li22_sym(Z3, pv, Z3, cv, mv))
DTv = {E: [x + y for x, y in zip(
        [Fr(vv) * ({-1:1,0:-1,1:0}[E]) for vv in to81(S22m)],
        [Fr(vv) * ({-1:-1,0:-1,1:0}[E]) for vv in to81(S22p)])] for E in (-1,0,1)}
Rv = {E: [x - y for x, y in zip(OMv[E], DTv[E])] for E in (-1,0,1)}

li4v = {s: to81(li_sym(4, rv, omv)) for s, rv, omv in nonconst}
logs = {'c': cv, '1-c': mv, '1+c': pv}
prodv = {}
for s, rv, omv in nonconst:
    S3 = li_sym(3, rv, omv); S2 = li_sym(2, rv, omv)
    for hn, hv in logs.items():
        prodv[f"Li3({s})*log({hn})"] = to81(sshuffle(S3, log_sym(hv)))
        for h2n, h2v in logs.items():
            prodv[f"Li2({s})*log({hn})*log({h2n})"] = to81(sshuffle(S2, sshuffle(log_sym(hv), log_sym(h2v))))
    prodv[f"Li2({s})**2"] = to81(sshuffle(S2, S2))
for i, (s1, rv1, omv1) in enumerate(nonconst):
    for s2, rv2, omv2 in nonconst[i+1:]:
        prodv[f"Li2({s1})*Li2({s2})"] = to81(sshuffle(li_sym(2, rv1, omv1), li_sym(2, rv2, omv2)))
for combo in itertools.combinations_with_replacement(list(logs), 4):
    S = log_sym(logs[combo[0]])
    for h in combo[1:]:
        S = sshuffle(S, log_sym(logs[h]))
    prodv["*".join(f"log({h})" for h in combo)] = to81(S)
print("li4:", len(li4v), "prods:", len(prodv))
pickle.dump({'li4v': li4v, 'prodv': prodv, 'Rv': Rv, 'OMv': OMv, 'DTv': DTv}, open('solve81_v2.pkl','wb'))
# rank check
import numpy as np
p = (1<<31)-1
def rank_modp(vecs):
    basis = []
    for v in vecs:
        v = np.array([x % p for x in v], dtype=np.int64)
        for b, pw in basis:
            if v[pw]: v = (v - int(v[pw])*b) % p
        nz = np.nonzero(v)[0]
        if len(nz):
            pw = int(nz[0]); basis.append(((v*pow(int(v[pw]), p-2, p)) % p, pw))
    return len(basis)
allv = list(li4v.values()) + list(prodv.values())
ra = rank_modp(allv)
print("rank all:", ra, " R_E in span:", all(rank_modp(allv + [Rv[E]]) == ra for E in (-1,0,1)))
