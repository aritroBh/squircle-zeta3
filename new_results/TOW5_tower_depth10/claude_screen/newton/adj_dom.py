import sys
sys.dont_write_bytecode=True
W="<project>/research/08_theory/TOW5_tower_depth10"
sys.path.insert(0,W)
import sympy as sp
import prove_adj910_fast as P
p=P.p
for delta,s in [(9,3),(9,4),(10,4),(10,5)]:
    terms=[sp.cancel(P.eR(delta,s,P.d,p)+P.eR(delta,s+1,P.d,p))]+P.feed_terms(delta,s)
    bad=[];nd=0
    for t in terms:
        n_,d_=sp.fraction(sp.cancel(t))
        for f,e in sp.factor_list(d_)[1]:
            if f.has(p):
                nd+=1
                for r in sp.Poly(f,p).nroots(n=30):
                    if abs(sp.im(r))<1e-20 and sp.re(r)>=4-1e-9: bad.append((f,r))
    print((delta,s),"denominator factors seen",nd,"real roots >=4:",bad)
