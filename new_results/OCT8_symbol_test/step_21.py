import pickle, sympy as sp, itertools
from fractions import Fraction as Fr
import sys; sys.path.insert(0,".")
from oct8 import *

sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational
letters = [2,3,5]

# ---- named weight-2 tensors ----
cands = list(sols) + [Fr(2),Fr(3),Fr(5),Fr(1,2),Fr(1,3),Fr(1,5),Fr(-1,2),Fr(-1,3),Fr(-1,5),Fr(-2),Fr(-4)]
W2_named = []   # (kind, arg) kind in ('Li2', z) or ('sym', a, b)
for z in cands:
    if z in (0,1): continue
    W2_named.append(('Li2', z))
for a in letters:
    for b in letters:
        W2_named.append(('sym', a, b))
def W2_tensor(item):
    if item[0] == 'Li2':
        z = item[1]
        return tensor2(vecfun(1-z), vecfun(z))
    _, a, b = item
    return sadd2(tensor2({a:Fr(1)},{b:Fr(1)}), tensor2({b:Fr(1)},{a:Fr(1)}))
def tensor2(u, v):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            S[(a,b)] = S.get((a,b), Fr(0)) + ca*cb
    return {k: x for k,x in S.items() if x != 0}
def sadd2(*Ss): return sadd(*Ss)

bidx2 = {t:i for i,t in enumerate(itertools.product(letters, repeat=2))}
def coord2(T):
    w = sp.zeros(1, 9)
    for t, c in T.items():
        w[bidx2[t]] += sp.Rational(c)
    return w
M = sp.Matrix.hstack(*[coord2(W2_tensor(it)).T for it in W2_named])
piv = M.rref()[1]
W2_basis = [W2_named[j] for j in piv]
print("W2 named basis:", W2_basis)

# ---- Decomp named generators: a sh W2b ----
# function: Li2-basis elt (1-w)(x)w sh a  <-> -Li2(w)*log(a) ;  sym(b,c) sh a <-> log a log b log c
Q = Quotient(letters)
D_named = []  # (description, symbol)
for a in letters:
    for it in W2_basis:
        T = W2_tensor(it)  # dict (b,c)->coeff
        gensym = {}
        for (b,c), cf in T.items():
            for t in [(a,b,c),(b,a,c),(b,c,a)]:
                gensym[t] = gensym.get(t, Fr(0)) + cf
        if it[0]=='Li2':
            desc = ('prod', -Li2 := None)  # placeholder
            desc = f"-Li2({it[1]})*log{a}"
        else:
            desc = f"log{a}*log{it[1]}*log{it[2]}"
        D_named.append((desc, {k:v for k,v in gensym.items() if v != 0}))
print("# Decomp named:", len(D_named))

# ---- Li3 generators ----
G_named = [(f"Li3({z})", S_Li3(z, vecfun)) for z in sols]
gens = [Q.coord(S) for _, S in G_named]
Dn   = [Q.coord(S) for _, S in D_named]
bigM = sp.Matrix.hstack(*[g.T for g in gens] + [d.T for d in Dn])
print("rank [G|Dnamed] =", bigM.rank(), " (should be 27)")

# ---- targets ----
targets = {
 "Lam1": (Fr(-1,3), Fr(1,3)),
 "Lam2": (Fr(-1,9), Fr(3,1)),
 "Lam3": (Fr(1,4), Fr(2,5)),
 "Lam4": (Fr(1,10), Fr(5,2)),
}
TS = {n: S_Li21(x,y,vecfun) for n,(x,y) in targets.items()}
# first check: is target in Decomp alone? (Pi=0?)
DM = sp.Matrix.hstack(*[d.T for d in Dn])
for n, S in TS.items():
    t = Q.coord(S)
    print(n, "Pi S = 0 (in Decomp):", sp.Matrix.hstack(DM, t.T).rank() == DM.rank(),
          " rank[t] nonzero:", any(x != 0 for x in t))
sol_out = {}
for n, S in TS.items():
    t = Q.coord(S)
    c = bigM.gauss_jordan_solve(t.T)[0]
    c = [sp.Rational(x) for x in c]  # particular solution (taus->0)
    sol_out[n] = c
    nzG = [(G_named[i][0], c[i]) for i in range(len(G_named)) if c[i] != 0]
    nzD = [(D_named[i][0], c[len(G_named)+i]) for i in range(len(D_named)) if c[len(G_named)+i] != 0]
    print("==", n, ": Li3 terms:", len(nzG), " product terms:", len(nzD))
    print("   Li3:", nzG[:12])
pickle.dump({"sol_out":{n:[str(x) for x in c] for n,c in sol_out.items()},
             "G_names":[g for g,_ in G_named], "D_names":[d for d,_ in D_named]},
            open("symbol_identities.pkl","wb"))
print("saved symbol_identities.pkl")
