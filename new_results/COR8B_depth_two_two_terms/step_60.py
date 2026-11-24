
import sys, pickle, itertools
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import sympy as sp
br = pickle.load(open('basis_red81.pkl','rb'))['basis_red']
piv_rows = {}
for name, row in br:
    piv_rows[min(row)] = (row, name)
def resid_span(vec81):
    v = {i: x for i,x in enumerate(vec81) if x}
    while v:
        i = min(v)
        if i not in piv_rows: return v
        row, name = piv_rows[i]
        f = v[i]/row[i]
        for j,x in row.items():
            if j==i: continue
            nv = v.get(j, Fr(0)) - f*x
            if nv: v[j] = nv
            elif j in v: del v[j]
        del v[i]
    return v
c = sp.symbols('c')
def factor_vec(expr):
    expr = sp.simplify(expr)
    num, den = sp.fraction(sp.cancel(expr))
    vec = [0,0,0,0]
    for poly, sgn in [(num,1),(den,-1)]:
        poly = sp.Poly(poly, c)
        content, factors = sp.factor_list(poly)
        cont = sp.Rational(content);  cont = -cont if cont<0 else cont
        for prime, e in sp.factorint(int(cont.p)).items():
            if prime != 2: return None
            vec[0] += sgn*e
        for prime, e in sp.factorint(int(cont.q)).items():
            if prime != 2: return None
            vec[0] -= sgn*e
        for f, e in factors:
            f = sp.Poly(f, c)
            if f.degree() == 0: continue
            fc = f.as_expr(); matched = False
            for basis_expr, idx in [(c,1),(1-c,2),(1+c,3)]:
                ratio = sp.simplify(fc/basis_expr)
                if ratio.is_number:
                    r = sp.Rational(ratio); r = -r if r<0 else r
                    for prime, ee in sp.factorint(int(r.p)).items():
                        if prime != 2: return None
                        vec[0] += sgn*e*ee
                    for prime, ee in sp.factorint(int(r.q)).items():
                        if prime != 2: return None
                        vec[0] -= sgn*e*ee
                    vec[idx] += sgn*e; matched = True; break
            if not matched: return None
    return tuple(vec)
E = [(1,0,0),(0,1,0),(0,0,1)]
def expand(sym):
    out = {}
    for k,c0 in sym.items():
        opts = [[(i,vi) for i,vi in enumerate(letter) if vi] for letter in k]
        for combo in itertools.product(*opts):
            mult = Fr(c0)
            for i,vi in combo: mult *= vi
            nk = tuple(i for i,vi in combo)
            out[nk] = out.get(nk, Fr(0)) + mult
        out = {k:v for k,v in out.items() if v}
    return out
def vec81(sym):
    e = expand(sym)
    v = [Fr(0)]*81
    for k,c0 in e.items():
        v[k[0]*27 + k[1]*9 + k[2]*3 + k[3]] = c0
    return v
def shuffle(symA, symB):
    out = {}
    def rec(pref, A, B, cur):
        if not A and not B: out[pref] = out.get(pref, Fr(0)) + cur; return
        if A: rec(pref+(A[0],), A[1:], B, cur)
        if B: rec(pref+(B[0],), A, B[1:], cur)
    for ka, ca in symA.items():
        for kb, cb in symB.items():
            rec((), list(ka), list(kb), ca*cb)
    return {k:v for k,v in out.items() if v}
def S_Li_n(g, n):
    vg = factor_vec(g); v1g = factor_vec(1-g)
    u, x = v1g[1:], vg[1:]
    return {tuple([u]+[x]*(n-1)): Fr(-1)}
def S_log(g):
    vg = factor_vec(g)
    return {(vg[1:],): Fr(1)}   # S(log g) = [g]

# test COR8's used products
prods = [
    ('Li2((1-c)*(c+1))*log(c)*log(1-c)', ('Li2','(1-c)*(c+1)'), ('log','c'), ('log','1-c')),
    ('Li2((1-c)^2/(c+1)^2)*log(c)*log(c)', ('Li2','(1-c)**2/(c+1)**2'), ('log','c'), ('log','c')),
    ('Li2((1-c)/(c+1))*log(1-c)*log(1-c)', ('Li2','(1-c)/(c+1)'), ('log','1-c'), ('log','1-c')),
    ('Li3(-2c/(1-c))*log(1+c)', ('Li3','-2*c/(1-c)'), ('log','1+c')),
    ('Li3(-c)*log(1-c)', ('Li3','-c'), ('log','1-c')),
    ('Li2(-c)**2', ('Li2','-c'), ('Li2','-c')),
    ('Li3((1-c)/(c+1))*log(c)', ('Li3','(1-c)/(c+1)'), ('log','c')),
]
for name, *parts in prods:
    s = None
    for kind, a in parts:
        g = sp.sympify(a)
        if kind=='Li2': sp_s = S_Li_n(g,2)
        elif kind=='Li3': sp_s = S_Li_n(g,3)
        elif kind=='log': sp_s = S_log(g)
        s = sp_s if s is None else shuffle(s, sp_s)
    r = resid_span(vec81(s))
    print(("OK " if not r else "OUT"), name, ("" if not r else f"resid {len(r)} coords e.g. {sorted(r.items())[:4]}"))
