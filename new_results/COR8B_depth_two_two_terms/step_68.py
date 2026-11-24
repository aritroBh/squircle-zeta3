
import sys, pickle, itertools, re
sys.path.insert(0,'.')
from fractions import Fraction as Fr
import sympy as sp
c = sp.symbols('c')
exec(open('step_66_helpers.py').read()) if False else None
# rebuild minimal helpers (same as previous cell, but only needed defs)
def split_factors(s):
    parts, depth, cur = [], 0, ''
    for ch in s:
        if ch == '(': depth += 1; cur += ch
        elif ch == ')': depth -= 1; cur += ch
        elif ch == '*' and depth == 0:
            if cur: parts.append(cur); cur = ''
        else: cur += ch
    if cur: parts.append(cur)
    return parts
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
    return {tuple([v1g[1:]]+[vg[1:]]*(n-1)): Fr(-1)}
def S_log(g):
    return {(factor_vec(g)[1:],): Fr(1)}
def symbol_of_item(item):
    if item[0] == 'Li4':
        return S_Li_n(sp.sympify(item[1]), 4)
    s = item[1].strip()
    dup_last = False
    if s.endswith('**2'):
        s = s[:-3].strip(); dup_last = True
    parts = [fac.strip() for fac in split_factors(s)]
    if dup_last: parts.append(parts[-1])
    sym = None
    for fac in parts:
        m = re.match(r'(Li2|Li3|log)\((.*)\)$', fac)
        fname, arg = m.group(1), sp.sympify(m.group(2))
        fs_ = S_Li_n(arg, 2) if fname=='Li2' else S_Li_n(arg,3) if fname=='Li3' else S_log(arg)
        sym = fs_ if sym is None else shuffle(sym, fs_)
    return sym

fs = pickle.load(open('../COR8_depth_two_any_angle/full_sols_v3.pkl','rb'))
alln, sols = fs['alln'], fs['sols']
T = pickle.load(open('targets81.pkl','rb'))
cache = {}
for e in (-1,0):
    acc = [Fr(0)]*81
    for i, item in enumerate(alln):
        cf = sols[e][i]
        if not cf: continue
        key = str(item)
        if key not in cache: cache[key] = vec81(symbol_of_item(item))
        v = cache[key]
        for j in range(81): acc[j] += cf*v[j]
    diff = [(j, acc[j]-T[e][j]) for j in range(81) if acc[j]-T[e][j]]
    print(f"e={e} differing components (COR8sol - myTarget):")
    for j, d in diff:
        print(f"   idx {j} letters ({j//27},{(j//9)%3},{(j//3)%3},{j%3}): {d}")
