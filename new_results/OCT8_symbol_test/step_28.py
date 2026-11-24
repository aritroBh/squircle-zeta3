import pickle, sympy as sp, itertools
from fractions import Fraction as Fr
import sys; sys.path.insert(0,".")
from oct8 import *

sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational
letters = [2,3,5]

def tensor2(u, v):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            S[(a,b)] = S.get((a,b), Fr(0)) + ca*cb
    return {k: x for k,x in S.items() if x != 0}
cands = list(sols) + [Fr(2),Fr(3),Fr(5),Fr(1,2),Fr(1,3),Fr(1,5),Fr(-1,2),Fr(-1,3),Fr(-1,5),Fr(-2),Fr(-4)]
W2_named = [('Li2', z) for z in cands if z not in (0,1)] + [('sym', a, b) for a in letters for b in letters]
def W2_tensor(item):
    if item[0] == 'Li2':
        z = item[1]; return tensor2(vecfun(1-z), vecfun(z))
    _, a, b = item
    return sadd(tensor2({a:Fr(1)},{b:Fr(1)}), tensor2({b:Fr(1)},{a:Fr(1)}))
bidx2 = {t:i for i,t in enumerate(itertools.product(letters, repeat=2))}
M = sp.Matrix.hstack(*[ (lambda T: (lambda w: w)(sp.zeros(1,9)))(W2_tensor(it)).T for it in []]) if False else None
def coord2(T):
    w = sp.zeros(1, 9)
    for t, c in T.items(): w[bidx2[t]] += sp.Rational(c)
    return w
M = sp.Matrix.hstack(*[coord2(W2_tensor(it)).T for it in W2_named])
W2_basis = [W2_named[j] for j in M.rref()[1]]

Q = Quotient(letters)
D_named = []
for a in letters:
    for it in W2_basis:
        T = W2_tensor(it); gensym = {}
        for (b,c), cf in T.items():
            for t in [(a,b,c),(b,a,c),(b,c,a)]:
                gensym[t] = gensym.get(t, Fr(0)) + cf
        gensym = {k:v for k,v in gensym.items() if v != 0}
        desc = f"-Li2({it[1]})*log{a}" if it[0]=='Li2' else f"log{a}*log{it[1]}*log{it[2]}"
        D_named.append((desc, gensym))
G_named = [(f"Li3({z})", S_Li3(z, vecfun)) for z in sols]
gens = [Q.coord(S) for _, S in G_named]
Dn   = [Q.coord(S) for _, S in D_named]
bigM = sp.Matrix.hstack(*[g.T for g in gens] + [d.T for d in Dn])

def S_Li12_xy(x, y, vecfun):
    # S(Li_{1,2}(x,y)) = (xy-1)(x)[x(y-1)/(x-1)](x)y + (y-1)(x)(x-1)(x)y
    #                  + (xy-1)(x)(xy)(x)[x/(x-1)] + (y-1)(x)y(x)(x-1)
    L1 = vecfun(x*y-1); L2 = vecfun(x*(y-1)/(x-1)); L3 = vecfun(y)
    L4 = vecfun(y-1);   L5 = vecfun(x-1);           L6 = vecfun(x*y)
    L7 = vecfun(x/(x-1))
    return sadd(tensor3(L1,L2,L3), tensor3(L4,L5,L3), tensor3(L1,L6,L7), tensor3(L4,L3,L5))

# OCT7 convention-B constants: LamB(a,c) = standard Li_{1,2}(c, a)
targetsB = {
 "Lam1": (Fr(1,3), Fr(-1,3)),    # Li_{1,2}(1/3, -1/3)
 "Lam2": (Fr(3,1), Fr(-1,9)),    # Li_{1,2}(3, -1/9)
 "Lam3": (Fr(2,5), Fr(1,4)),     # Li_{1,2}(2/5, 1/4)
 "Lam4": (Fr(5,2), Fr(1,10)),    # Li_{1,2}(5/2, 1/10)
}
TSB = {n: S_Li12_xy(x,y,vecfun) for n,(x,y) in targetsB.items()}
# verify letters all in {2,3,5}
allL = set()
for S in TSB.values():
    for t in S: allL.update(t)
print("letters conv-B targets:", sorted(allL))

solB = {}
for n, S in TSB.items():
    t = Q.coord(S)
    c = bigM.gauss_jordan_solve(t.T)[0]
    c = [sp.Rational(x.subs({s:0 for s in x.free_symbols})) for x in c]
    solB[n] = c
    nzG = [(G_named[i][0], c[i]) for i in range(len(G_named)) if c[i] != 0]
    nzD = [(D_named[i][0], c[len(G_named)+i]) for i in range(len(D_named)) if c[len(G_named)+i] != 0]
    print("==", n, ": #Li3:", len(nzG), " #prod:", len(nzD))
    print("   Li3:", nzG)
    print("   prd:", nzD)
pickle.dump({"sol_out":{n:[str(x) for x in c] for n,c in solB.items()},
             "G_names":[g for g,_ in G_named], "D_names":[d for d,_ in D_named]},
            open("symbol_identities_B.pkl","wb"))
print("saved symbol_identities_B.pkl")
