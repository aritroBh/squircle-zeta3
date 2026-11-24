
import sys, pickle, itertools
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import sympy as sp
br = pickle.load(open('basis_red81.pkl','rb'))['basis_red']
piv_rows = {}
for name, row in br:
    piv_rows[min(row)] = (row, name)
def in_span(vec81):
    v = {i: x for i,x in enumerate(vec81) if x}
    while v:
        i = min(v)
        if i not in piv_rows: return False
        row, name = piv_rows[i]
        f = v[i]/row[i]
        for j,x in row.items():
            if j==i: continue
            nv = v.get(j, Fr(0)) - f*x
            if nv: v[j] = nv
            elif j in v: del v[j]
        del v[i]
    return True

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

# COR8 Li4 args from the e=-1 used list
args = ['(1-c)*(c+1)','(1-c)**2/(c+1)**2','(c+1)**2/(1-c)**2','-2*c/(1-c)','-4*c/(1-c)**2','-c',
        '-c**2/((1-c)*(c+1))','-c/(1-c)','1-c','1/((1-c)*(c+1))','1/(1-c)','1/(c+1)','1/2-c/2',
        '2*c/(c+1)','2/(c+1)','4*c/(c+1)**2','c','c+1','c/(c+1)','(1-c)/(c+1)','(c-1)/(c+1)','-c/(c+1)']
for a in args:
    g = sp.sympify(a)
    vg = factor_vec(g); v1g = factor_vec(1-g)
    if vg is None or v1g is None:
        print(f"{a}: letters OUTSIDE alphabet")
        continue
    u = v1g[1:]; x = vg[1:]
    s = {(u,x,x,x): Fr(-1)}
    ok = in_span(vec81(s))
    if not ok: print(f"Li4({a}): NOT in span   ([1-g]={u}, [g]={x})")
print("done checking Li4 args")
