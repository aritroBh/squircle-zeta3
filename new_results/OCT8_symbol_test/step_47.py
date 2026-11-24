import mpmath as mp, pickle, sympy as sp, itertools, re, sys
from fractions import Fraction as Fr
sys.path.insert(0,".")
from oct8 import *
mp.mp.dps = 300
def M(q):
    q = Fr(q); return mp.mpf(q.numerator)/q.denominator
sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational; letters = [2,3,5]
Q = Quotient(letters)
def tensor2(u, v):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            S[(a,b)] = S.get((a,b), Fr(0)) + ca*cb
    return {k: x for k,x in S.items() if x != 0}
def W2_tensor(item):
    if item[0]=='Li2': return tensor2(vecfun(1-item[1]), vecfun(item[1]))
    _,a,b = item; return sadd(tensor2({a:Fr(1)},{b:Fr(1)}), tensor2({b:Fr(1)},{a:Fr(1)}))
cands = list(sols) + [Fr(2),Fr(3),Fr(5),Fr(1,2),Fr(1,3),Fr(1,5),Fr(-1,2),Fr(-1,3),Fr(-1,5),Fr(-2),Fr(-4)]
W2_named = [('Li2', z) for z in cands if z not in (0,1)] + [('sym', a, b) for a in letters for b in letters]
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
G_list = [(z, S_Li3(z, vecfun)) for z in sols if z != Fr(-3)]  # EXCLUDE Li3(-3)
gens = [Q.coord(S) for _,S in G_list]; Dn = [Q.coord(S) for _,S in D_named]
bigM = sp.Matrix.hstack(*[g.T for g in gens]+[d.T for d in Dn])
targ = S_Li3(Fr(-3), vecfun)
t = Q.coord(targ)
c = bigM.gauss_jordan_solve(t.T)[0]
c = [sp.Rational(x.subs({s:0 for s in x.free_symbols})) for x in c]
logv = {p: mp.log(p) for p in (2,3,5)}
tot = mp.mpf(0)
for i,(z,_) in enumerate(G_list):
    if c[i]: tot += M(c[i])*mp.polylog(3, M(z))
for j,(dn,_) in enumerate(D_named):
    cf = c[len(G_list)+j]
    if not cf: continue
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn); w,p = Fr(m.group(1)), int(m.group(2))
    tot += M(cf)*(-mp.polylog(2, M(w)))*logv[p]
F = mp.polylog(3, -3)
rem = F - tot
print("Li3(-3) =", mp.nstr(F, 15), " recovered =", mp.nstr(tot,15), " rem =", mp.nstr(abs(rem),3))
rel = mp.pslq([rem, mp.zeta(3), mp.zeta(2)*mp.log(2), mp.zeta(2)*mp.log(3), mp.zeta(2)*mp.log(5)],
              tol=mp.mpf('1e-260'), maxcoeff=10**60, maxsteps=1000000)
print("PSLQ:", rel)
print("nonzero coeffs used:", [(G_list[i][0], c[i]) for i in range(len(G_list)) if c[i]][:8], "...",
      [(D_named[j][0], c[len(G_list)+j]) for j in range(len(D_named)) if c[len(G_list)+j]])
