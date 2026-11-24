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
        z = item[1]
        return tensor2(vecfun(1-z), vecfun(z))
    _, a, b = item
    return sadd(tensor2({a:Fr(1)},{b:Fr(1)}), tensor2({b:Fr(1)},{a:Fr(1)}))

bidx2 = {t:i for i,t in enumerate(itertools.product(letters, repeat=2))}
def coord2(T):
    w = sp.zeros(1, 9)
    for t, c in T.items():
        w[bidx2[t]] += sp.Rational(c)
    return w
M = sp.Matrix.hstack(*[coord2(W2_tensor(it)).T for it in W2_named])
piv = M.rref()[1]
W2_basis = [W2_named[j] for j in piv]

Q = Quotient(letters)
D_named = []
for a in letters:
    for it in W2_basis:
        T = W2_tensor(it)
        gensym = {}
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

targets = {
 "Lam1": (Fr(-1,3), Fr(1,3)),
 "Lam2": (Fr(-1,9), Fr(3,1)),
 "Lam3": (Fr(1,4), Fr(2,5)),
 "Lam4": (Fr(1,10), Fr(5,2)),
}
TS = {n: S_Li21(x,y,vecfun) for n,(x,y) in targets.items()}

sol_out = {}
for n, S in TS.items():
    t = Q.coord(S)
    c = bigM.gauss_jordan_solve(t.T)[0]
    c = [sp.Rational(x.subs({s:0 for s in x.free_symbols})) for x in c]
    sol_out[n] = c
    nzG = [(G_named[i][0], c[i]) for i in range(len(G_named)) if c[i] != 0]
    nzD = [(D_named[i][0], c[len(G_named)+i]) for i in range(len(D_named)) if c[len(G_named)+i] != 0]
    print("==", n, ": #Li3 terms:", len(nzG), " #product terms:", len(nzD))
    print("   Li3:", nzG)
    print("   prd:", nzD[:10])
pickle.dump({"sol_out":{n:[str(x) for x in c] for n,c in sol_out.items()},
             "G_names":[g for g,_ in G_named], "D_names":[d for d,_ in D_named]},
            open("symbol_identities.pkl","wb"))
print("saved")
with open("progress.log","a") as f:
    f.write("STEP Q-span: Li3(99 S-sols)+named Decomp spans all of V^3 (rank 27). All four Lams IN span. "
            "Explicit identities extracted to symbol_identities.pkl (particular sol, taus=0).\n")
