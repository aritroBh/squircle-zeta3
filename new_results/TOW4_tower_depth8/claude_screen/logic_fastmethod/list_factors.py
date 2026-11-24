import sys, sympy as sp
import prove_adj78_fast as A
p=A.p
delta,s=int(sys.argv[1]),int(sys.argv[2])
lhs = sp.cancel(A.eR(delta, s, A.d, p) + A.eR(delta, s+1, A.d, p))
terms = A.feed_terms(delta, s)
exprs=[lhs]+[-t for t in terms]
for j,e in enumerate(exprs):
    n,dn=sp.fraction(sp.cancel(e))
    if n==0: continue
    c,fl=sp.factor_list(dn)
    print(j, 'content',c, [(str(sp.factor(f)),ee) for f,ee in fl])
