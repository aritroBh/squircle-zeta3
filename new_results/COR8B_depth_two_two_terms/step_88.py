
import sys, pickle, itertools
sys.path.insert(0,'.')
from cor8b_sym import *
from fractions import Fraction as Fr
import sympy as sp
c = sp.symbols('c')

# HPL weight-4 words at S-unit args: letters [g],[1-g],[1+g] must be in alphabet
vals = pickle.load(open('sunit_args_fixed.pkl','rb'))
def fv_plus(expr):
    # factor 1+g
    pass
gmapf = {0: 'x', 1: 'u', -1: 'w'}   # letters: [g],[1-g],[1+g]
words = list(itertools.product([0,1,-1], repeat=4))
dataH = []
used_args = 0
for arg_s, (vt, v1t) in vals.items():
    g = sp.sympify(arg_s)
    # need [1+g] in alphabet too
    import sympy as sp_
    # factor 1+g
    num, den = sp.fraction(sp.cancel(1+g))
    ok = True
    vec = [0,0,0,0]
    for poly, sgn in [(num,1),(den,-1)]:
        poly = sp.Poly(poly, c)
        content, factors = sp.factor_list(poly)
        cont = sp.Rational(content);  cont = -cont if cont<0 else cont
        for prime, e in sp.factorint(int(cont.p)).items():
            if prime != 2: ok = False; break
            vec[0] += sgn*e
        for prime, e in sp.factorint(int(cont.q)).items():
            if prime != 2: ok = False; break
            vec[0] -= sgn*e
        if not ok: break
        for f, e in factors:
            f = sp.Poly(f, c)
            if f.degree() == 0: continue
            fc = f.as_expr(); matched = False
            for basis_expr, idx in [(c,1),(1-c,2),(1+c,3)]:
                ratio = sp.simplify(fc/basis_expr)
                if ratio.is_number:
                    r = sp.Rational(ratio); r = -r if r<0 else r
                    for prime, ee in sp.factorint(int(r.p)).items():
                        if prime != 2: ok=False; break
                        vec[0] += sgn*e*ee
                    for prime, ee in sp.factorint(int(r.q)).items():
                        if prime != 2: ok=False; break
                        vec[0] -= sgn*e*ee
                    vec[idx] += sgn*e; matched = True; break
            if not matched: ok = False; break
        if not ok: break
    if not ok: continue
    used_args += 1
    gm = {0: vt, 1: v1t, -1: tuple(vec)}
    for w in words:
        sym = {tuple(gm[a] for a in reversed(w)): Fr(1)}
        dd = delta22(sym)
        if dd:
            dataH.append((str(g), ''.join(map(str,w)), dd))
print("args with [1+g] in alphabet:", used_args, " HPL deltas:", len(dataH))
pickle.dump(dataH, open('hpl_delta_all.pkl','wb'))
# combine with MPL database and check span
dataM = pickle.load(open('pairs_delta_all_v2.pkl','rb'))
USE = [(1,3),(1,4),(2,3),(2,4),(3,4),(3,5),(4,5)]
ALLCOMP = sorted({(i,j) for i in range(6) for j in range(6) if i<j})
cidx = {k:i for i,k in enumerate(ALLCOMP)}
def vec_of(dd):
    v = [Fr(0)]*len(ALLCOMP)
    for k,val in dd.items():
        kk=(min(k),max(k)); sgn=1 if k==kk else -1
        v[cidx[kk]] += sgn*val
    return v
vecs = []
for (x0,y0,d22,d31,d13) in dataM:
    for dd in (d22,d31,d13):
        vecs.append(tuple(vec_of(dd)))
for (g,w,dd) in dataH:
    vecs.append(tuple(vec_of(dd)))
def rank(rows):
    rows = [list(r) for r in rows if any(r)]
    if not rows: return 0
    m,n=len(rows),len(rows[0]); r=0
    for col in range(n):
        piv=None
        for i in range(r,m):
            if rows[i][col]: piv=i;break
        if piv is None: continue
        rows[r],rows[piv]=rows[piv],rows[r]
        pv=rows[r][col]
        for i in range(m):
            if i!=r and rows[i][col]:
                f=rows[i][col]/pv
                rows[i]=[a-f*b for a,b in zip(rows[i],rows[r])]
        r+=1
        if r==m: break
    return r
T0 = [Fr(0)]*len(ALLCOMP); T1 = [Fr(0)]*len(ALLCOMP)
T0[cidx[(3,5)]] = 1; T0[cidx[(4,5)]] = -1
T1[cidx[(3,5)]] = -1; T1[cidx[(4,5)]] = -1
R = rank(vecs)
print("combined span rank:", R)
print("T0 in span:", rank(vecs+[tuple(T0)])==R, " T1 in span:", rank(vecs+[tuple(T1)])==R)
