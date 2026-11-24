import sys, time
sys.dont_write_bytecode=True
W="<project>/research/08_theory/TOW5_tower_depth10"
sys.path.insert(0,W)
import sympy as sp
import prove_E3b910 as E
p=E.p
delta=int(sys.argv[1]); t0=time.time()
terms=E.feed_terms(delta)
ok=E.check_zero(terms)
print(f"E3 delta={delta}: exact numerator==0: {ok} ({time.time()-t0:.1f}s)")
# domain: integer roots >=4 (and any real root >= 4) of every cancelled denominator used
dens=set()
for cache in (E._ehat,E._D,E._Dc,E._g):
    for k,v in cache.items():
        n_,d_=sp.fraction(sp.cancel(v)); dens.add(sp.factor(d_))
for t in terms:
    n_,d_=sp.fraction(sp.cancel(t)); dens.add(sp.factor(d_))
bad=[]
for dn in dens:
    for f,e in sp.factor_list(dn)[1]:
        if f.has(p):
            for r in sp.Poly(f,p).nroots(n=30) if sp.degree(f,p)>0 else []:
                if abs(sp.im(r))<1e-20 and sp.re(r)>=4-1e-9: bad.append((f,r))
print("distinct denominators:",len(dens)); print("denominator roots with real part >= 4:",bad)
print("sample denominators:",list(dens)[:4])
