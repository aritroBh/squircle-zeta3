import mpmath as mp, pickle, re, sympy as sp, itertools, sys
from fractions import Fraction as Fr
sys.path.insert(0,".")
from oct8 import *
mp.mp.dps = 100
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator

sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational; letters = [2,3,5]
Q = Quotient(letters)

# Synthetic: F = 2*Li3(-3) + Li3(-5) - Li2(-3)*log2   (weight-3, reduces by construction)
# S(F) = 2*S(Li3(-3)) + S(Li3(-5)) - S(Li2(-3)log2)
SF = sadd(smul(S_Li3(Fr(-3),vecfun),2), S_Li3(Fr(-5),vecfun),
          smul(tensor3(vecfun(1-Fr(-3)), vecfun(-3), vecfun(2)), 1))  # -Li2(-3)*log2 symbol = +((1-w)(x)w) sh log2? 
# careful: S(-Li2(w) log a) = ((1-w)(x)w) sh a. Let me build it as shuffle:
def shuffle_a_w(a, w):
    T = tensor2(vecfun(1-w), vecfun(w))
    r = {}
    for (b,c),cf in T.items():
        for t in [(a,b,c),(b,a,c),(b,c,a)]:
            r[t] = r.get(t,Fr(0))+cf
    return r
SF = sadd(smul(S_Li3(Fr(-3),vecfun),2), S_Li3(Fr(-5),vecfun), shuffle_a_w(2, Fr(-3)))
# numeric F
Fnum = 2*mp.polylog(3,-3) + mp.polylog(3,-5) - mp.polylog(2,-3)*mp.log(2)

# decompose S(F) in [G|D]
cands = list(sols) + [Fr(2),Fr(3),Fr(5),Fr(1,2),Fr(1,3),Fr(1,5),Fr(-1,2),Fr(-1,3),Fr(-1,5),Fr(-2),Fr(-4)]
W2_named = [('Li2', z) for z in cands if z not in (0,1)] + [('sym', a, b) for a in letters for b in letters]
def tensor2(u, v):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            S[(a,b)] = S.get((a,b), Fr(0)) + ca*cb
    return {k: x for k,x in S.items() if x != 0}
def W2_tensor(item):
    if item[0]=='Li2': return tensor2(vecfun(1-item[1]), vecfun(item[1]))
    _,a,b = item; return sadd(tensor2({a:Fr(1)},{b:Fr(1)}), tensor2({b:Fr(1)},{a:Fr(1)}))
bidx2 = {t:i for i,t in enumerate(itertools.product(letters, repeat=2))}
def coord2(T):
    w = sp.zeros(1,9)
    for t,c in T.items(): w[bidx2[t]] += sp.Rational(c)
    return w
M2 = sp.Matrix.hstack(*[coord2(W2_tensor(it)).T for it in W2_named])
W2_basis = [W2_named[j] for j in M2.rref()[1]]
D_named = []
for a in letters:
    for it in W2_basis:
        gensym = {}
        for (b,c),cf in W2_tensor(it).items():
            for t in [(a,b,c),(b,a,c),(b,c,a)]:
                gensym[t] = gensym.get(t,Fr(0))+cf
        gensym = {k:v for k,v in gensym.items() if v != 0}
        desc = f"-Li2({it[1]})*log{a}" if it[0]=='Li2' else f"log{a}*log{it[1]}*log{it[2]}"
        D_named.append((desc, gensym))
G_named = [(f"Li3({z})", S_Li3(z, vecfun)) for z in sols]
gens = [Q.coord(S) for _,S in G_named]; Dn = [Q.coord(S) for _,S in D_named]
bigM = sp.Matrix.hstack(*[g.T for g in gens]+[d.T for d in Dn])
t = Q.coord(SF)
c = bigM.gauss_jordan_solve(t.T)[0]
c = [sp.Rational(x.subs({s:0 for s in x.free_symbols})) for x in c]
logv = {p: mp.log(p) for p in (2,3,5)}
sG = sum(M(c[i])*mp.polylog(3, M(Fr(G_names if False else G_named[i][0][4:-1]))) for i in range(len(G_named)) if c[i])
sD = mp.mpf(0)
for j,dn in enumerate(D_named):
    cf = c[len(G_named)+j]
    if not cf: continue
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn[0]); w,p = Fr(m.group(1)), int(m.group(2))
    sD += M(cf)*(-mp.polylog(2, M(w)))*logv[p]
print("F num:", mp.nstr(Fnum, 20))
print("recovered Li3 sum:", mp.nstr(sG,20))
print("recovered prod sum:", mp.nstr(sD,20))
print("diff (F - recovered):", mp.nstr(abs(Fnum - sG - sD),3))
# what ARE the recovered coefficients?
print("nonzero G coeffs:", [(G_named[i][0], c[i]) for i in range(len(G_named)) if c[i]])
print("nonzero D coeffs:", [(D_named[j][0], c[len(G_named)+j]) for j in range(len(D_named)) if c[len(G_named)+j]])
