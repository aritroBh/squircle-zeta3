import mpmath as mp, pickle, re, sympy as sp, itertools, sys
from fractions import Fraction as Fr
sys.path.insert(0,".")
from oct8 import *
mp.mp.dps = 150
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
G_named = [(f"Li3({z})", S_Li3(z, vecfun)) for z in sols]
Dn = [Q.coord(S) for _,S in D_named]
DM = sp.Matrix.hstack(*[d.T for d in Dn])

# choose 8 Li3 forming a Q-basis: try small-height z's
zpick = [Fr(-2),Fr(-3),Fr(-4),Fr(-5),Fr(-8),Fr(-9),Fr(-15),Fr(-24)]
Gpick = [(f"Li3({z})", S_Li3(z, vecfun)) for z in zpick]
Mpick = sp.Matrix.hstack(*[Q.coord(S).T for _,S in Gpick] + [d.T for d in Dn])
print("rank [Gpick|D]:", Mpick.rank(), "(want 27)")
solB = {}
targetsB = {"Lam1": (Fr(1,3), Fr(-1,3)), "Lam2": (Fr(3,1), Fr(-1,9)),
            "Lam3": (Fr(2,5), Fr(1,4)),  "Lam4": (Fr(5,2), Fr(1,10))}
def S_Li12_xy(x, y):
    L1 = vecfun(x*y-1); L2 = vecfun(x*(y-1)/(x-1)); L3 = vecfun(y)
    L4 = vecfun(y-1);   L5 = vecfun(x-1);           L6 = vecfun(x*y)
    L7 = vecfun(x/(x-1))
    return sadd(tensor3(L1,L2,L3), tensor3(L4,L5,L3), tensor3(L1,L6,L7), tensor3(L4,L3,L5))
TSB = {n: S_Li12_xy(x,y) for n,(x,y) in targetsB.items()}
for n,S in TSB.items():
    t = Q.coord(S)
    c = Mpick.gauss_jordan_solve(t.T)[0]  # unique now (square)
    c = [sp.Rational(x.subs({s:0 for s in x.free_symbols})) for x in c]
    solB[n] = c
    print("==", n)
    print("   Li3:", [(Gpick[i][0], c[i]) for i in range(8) if c[i] != 0])
    print("   prd:", [(D_named[j][0], c[8+j]) for j in range(len(D_named)) if c[8+j] != 0])
pickle.dump({"sol":{n:[str(x) for x in c] for n,c in solB.items()},
             "G_names":[g for g,_ in Gpick], "D_names":[d for d,_ in D_named],
             "zpick":[str(z) for z in zpick]},
            open("symbol_identities_B8.pkl","wb"))
print("saved symbol_identities_B8.pkl")
