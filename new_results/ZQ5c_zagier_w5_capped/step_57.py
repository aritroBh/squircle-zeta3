import json, math
import sympy as sp
from fractions import Fraction as Ff
G = [[int(x) for x in a] for a in json.load(open('out/d13/genuine_lattice_cap.json'))['G']]
names = json.load(open('out/d13/genuine_lattice_cap.json'))['atom_names']
print("genuine lattice rank:", len(G))
Tt = [a for a in G if a[0]!=0]
T0 = [a for a in G if a[0]==0]
print("T relations:", len(Tt), "T-free:", len(T0))
print("T-coefs:", [a[0] for a in Tt])
g=36771482407
print("factor(g):", sp.factorint(g))
print("g/13^5 =", Ff(g, 13**5), "  g mod 13:", g%13)
# find min T-coef vector: solve integer combo of Tt rows to get a[0]=g, then reduce by T0
# Use sympy HNF on the matrix of T-coefs
import numpy as np
# Build: we want combo c (len len(G)) with (sum c_i G_i)[0] = g and small norm.
# First find any vector with a[0]=g via extended gcd on the T-coefs of all genuine rows
coefs=[a[0] for a in G]
# extended gcd
def xgcd_list(vals):
    # returns combo u with sum u_i*vals_i = gcd
    import math
    M=sp.Matrix([[sp.Integer(v)] for v in vals]).T  # 1 x n
    # use smith normal form approach: gcd
    g0=0; u=[0]*len(vals)
    for i,v in enumerate(vals):
        if v==0: continue
        if g0==0: g0=v; u[i]=1
        else:
            gg,s,t=sp.gcdex(g0,v)
            u=[s*x for x in u]; u[i]=t; g0=gg
    return g0,u
g0,u=xgcd_list(coefs)
print("xgcd gives gcd", g0)
w=[sum(u[i]*G[i][j] for i in range(len(G))) for j in range(len(G[0]))]
print("min-T vector w: a[0]=", w[0], "maxcoef", max(abs(x) for x in w), "nterms", sum(1 for x in w if x))
json.dump(dict(w=[str(x) for x in w]), open('out/d13/minT_vec.json','w'))
