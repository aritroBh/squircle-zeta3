import pickle, sympy as sp, re, sys
from fractions import Fraction as Fr
sys.path.insert(0,".")
from oct8 import *
sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational; letters = [2,3,5]
Q = Quotient(letters)
dB = pickle.load(open("symbol_identities_B.pkl","rb"))
G_names, D_names = dB["G_names"], dB["D_names"]
G_syms = [S_Li3(Fr(gn[4:-1]), vecfun) for gn in G_names]
def tensor2(u, v):
    S = {}
    for a, ca in u.items():
        for b, cb in v.items():
            S[(a,b)] = S.get((a,b), Fr(0)) + ca*cb
    return {k: x for k,x in S.items() if x != 0}
def Dsym(dn):
    m = re.match(r"-Li2\((.*)\)\*log(\d+)", dn)
    w, p = Fr(m.group(1)), int(m.group(2))
    T = tensor2(vecfun(1-w), vecfun(w)); a = p
    r = {}
    for (b,c_), cf in T.items():
        for t in [(a,b,c_),(b,a,c_),(b,c_,a)]:
            r[t] = r.get(t, Fr(0)) + cf
    return {k:v for k,v in r.items() if v != 0}
D_syms = [Dsym(dn) for dn in D_names]
gens = [Q.coord(S) for S in G_syms]; Dn = [Q.coord(S) for S in D_syms]
bigM = sp.Matrix.hstack(*[g.T for g in gens] + [d.T for d in Dn])
NS = bigM.nullspace()
print("nullity:", len(NS))
pickle.dump({"NS":[[str(x) for x in nv] for nv in NS]}, open("nullspace.pkl","wb"))
print("saved nullspace.pkl")
