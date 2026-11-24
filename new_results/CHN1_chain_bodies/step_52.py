import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
from fractions import Fraction as F
import pickle, sympy as sp, engine as E, mzv_exact as M
forms={k:{tuple(int(x) for x in e.split(',')):F(c) for e,c in v.items()} for k,v in pickle.load(open('forms.pkl','rb')).items()}
x,y,z,w,s,t=sp.symbols('x y z w s t')
def to_sp(p):
    return sum(sp.Rational(c.numerator,c.denominator)*x**e[0]*y**e[1]*z**e[2]*w**e[3] for e,c in p.items())
# value vector of c_W on Pi1: monomial -> poly(s,t)
def value_vec(W, sub):
    T=E.table(W); out={}
    for b,mon in enumerate(T.monomials):
        cf=sp.expand(to_sp(forms[(W,b)]).subs(sub))
        if cf!=0: out[mon]=cf
    return out
# Gamma product series GP1 = G(1+s)G(1+s+t)/G(1+2s+t), GP2 = swap s<->t
def gamma_series_2var(forms_list, Wmax):
    # forms_list: list of (eps, (alpha,beta)) meaning Gamma(1+(alpha s+ beta t)u)^eps ; returns dict W-> {monomial->poly}
    L={}
    for eps,(al,be) in forms_list:
        for m in range(2,Wmax+1):
            L[m]=L.get(m,0)+eps*sp.Rational((-1)**m,m)*(al*s+be*t)**m
    G=[{} for _ in range(Wmax+1)]; G[0]={():sp.Integer(1)}
    for W in range(2,Wmax+1):
        acc={}
        for i in range(2,W+1):
            for mon,c in G[W-i].items():
                for mon2,c2 in M.stuffle((i,),mon).items():
                    # value: multiply letter (i,) stuffle into monomial 'mon' -- but monomials are PRODUCTS; need word-stuffle:
                    pass
        G[W]=acc
    return None
# NOTE: stuffing letters into a product monomial requires the value-level stuffle:
# the gamma series coefficients are sums of products of single zetas: represent monomials directly by concatenation (no stuffle):
def gamma_series_prod(L, Wmax):
    # L: dict m->poly; series = exp(sum L[m] zeta(m) u^m), monomials by concatenation
    G=[{} for _ in range(Wmax+1)]; G[0]={():sp.Integer(1)}
    for W in range(2,Wmax+1):
        acc={}
        for i in range(2,W+1):
            for mon,c in G[W-i].items():
                m2=tuple(sorted(mon+((i,),)))
                acc[m2]=acc.get(m2,0)+sp.Rational(i,W)*L[i]*c
        G[W]={m:sp.expand(c) for m,c in acc.items() if c!=0}
    return G
def L_from_forms(forms_list, Wmax):
    L={}
    for eps,(al,be) in forms_list:
        for m in range(2,Wmax+1):
            L[m]=sp.expand(L.get(m,sp.Integer(0))+sp.Rational(eps*((-1)**m),m)*(al*s+be*t)**m)
    return L
GP1f=[(1,(1,0)),(1,(1,1)),(-1,(2,1))]
GP2f=[(1,(0,1)),(1,(1,1)),(-1,(1,2))]
GP1=gamma_series_prod(L_from_forms(GP1f,11),11)
GP2=gamma_series_prod(L_from_forms(GP2f,11),11)
sub={x:s,y:s,z:t,w:t}
cV={W:value_vec(W,sub) for W in range(2,12)}
# test c_W = R1*GP1_W + R2*GP2_W (+ R3 if W==0) with rational functions R1,R2
# collect equations: for each W, each monomial: c = R1 g1 + R2 g2
R1,R2=sp.symbols('R1 R2')
eqs=[]
for W in range(2,12):
    keys=set(cV[W])|set(GP1[W])|set(GP2[W])
    for mon in keys:
        eqs.append((W,mon, sp.expand(cV[W].get(mon,0) - R1*GP1[W].get(mon,0) - R2*GP2[W].get(mon,0))))
# solve from W=2,3 equations (R1,R2 rational funcs): take 2 equations
sols=sp.solve([eqs[0][2], eqs[1][2]],[R1,R2], dict=True)
print('candidate R1,R2 from first two eqs:', sols)
if sols:
    sol=sols[0]
    bad=0
    for W,mon,e in eqs:
        if sp.expand(e.subs(sol))!=0:
            bad+=1
            if bad<=5: print('  FAIL at W=%d mon=%s: %s'%(W,mon,sp.factor(e.subs(sol))))
    print('two-term GP1,GP2 test: %d of %d equations fail'%(bad,len(eqs)))
print('t=%.1fs'%(time.time()-t0))
