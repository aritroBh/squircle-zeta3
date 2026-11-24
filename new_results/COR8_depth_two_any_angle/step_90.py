
import pickle, sys, json, itertools
sys.path.insert(0, "<project>/research/08_theory/COR8_depth_two_any_angle")
from fractions import Fraction as Fr
import sympy as sp
from cor8_sym import *
c = sp.Symbol('c')
C4 = "<project>/research/08_theory/COR4_M3alpha_classical/"
Z3 = (Fr(0),)*3
def v3(j, val=1):
    v = list(Z3); v[j] = Fr(val); return tuple(v)
cv, mv, pv = v3(0), v3(1), v3(2)

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
sa = [(sp.sstr(r), v_r[1:], v_om[1:]) for (v_r, v_om), r in sorted(args.items(), key=lambda kv: str(kv[1]))]
polefree = [(s, rv, omv) for s, rv, omv in sa if rv != Z3 and sp.sympify(s).subs(c, 0).is_finite]
print(len(polefree), "pole-free nonconstant args")
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
li4v = {s: to81(li_sym(4, rv, omv)) for s, rv, omv in polefree}
logs = {'c': cv, '1-c': mv, '1+c': pv}
prodv = {}
for s, rv, omv in polefree:
    S3 = li_sym(3, rv, omv); S2 = li_sym(2, rv, omv)
    for hn, hv in logs.items():
        prodv[f"Li3({s})*log({hn})"] = to81(sshuffle(S3, log_sym(hv)))
        for h2n, h2v in logs.items():
            prodv[f"Li2({s})*log({hn})*log({h2n})"] = to81(sshuffle(S2, sshuffle(log_sym(hv), log_sym(h2v))))
    prodv[f"Li2({s})**2"] = to81(sshuffle(S2, S2))
for i, (s1, rv1, omv1) in enumerate(polefree):
    for s2, rv2, omv2 in polefree[i+1:]:
        prodv[f"Li2({s1})*Li2({s2})"] = to81(sshuffle(li_sym(2, rv1, omv1), li_sym(2, rv2, omv2)))
for combo in itertools.combinations_with_replacement(list(logs), 4):
    S = log_sym(logs[combo[0]])
    for h in combo[1:]:
        S = sshuffle(S, log_sym(logs[h]))
    prodv["*".join(f"log({h})" for h in combo)] = to81(S)
print("li4:", len(li4v), "prods:", len(prodv))
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
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/solve81_v2.pkl','rb'))
Rv = D['Rv']
allv = list(li4v.values()) + list(prodv.values())
ra = rank_modp(allv)
print("rank:", ra, "R_E in span:", all(rank_modp(allv + [Rv[E]]) == ra for E in (-1,0,1)))
pickle.dump({'li4v': li4v, 'prodv': prodv, 'Rv': Rv}, open('solve81_v3.pkl','wb'))
