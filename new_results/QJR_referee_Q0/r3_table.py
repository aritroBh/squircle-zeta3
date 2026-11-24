# Referee QJR: the 15 table instances, re-derived from scratch (own code).
from rcommon import *
import json
P2=pi**2
def N_(x): return float(sp.N(x,30))
def refl(x):   # Li2(x)+Li2(1-x)=pi^2/6-log x log(1-x), 0<x<1
    assert 0<N_(x)<1; return [(1,x),(1,1-x)], P2/6-logS(x)*logS(1-x)
def inv(x):    # Li2(x)+Li2(1/x)=-pi^2/6-1/2 log^2(-x), x<0
    assert N_(x)<0; return [(1,x),(1,1/x)], -P2/6-logS(-x)**2/2
def land(x):   # Li2(x)+Li2(x/(x-1))=-1/2 log^2(1-x), x<1
    assert N_(x)<1; return [(1,x),(1,x/(x-1))], -logS(1-x)**2/2
def dup(x):    # Li2(x)+Li2(-x)-1/2 Li2(x^2)=0, |x|<=1
    assert abs(N_(x))<=1; return [(1,x),(1,-x),(-sp.Rational(1,2),x**2)], 0
def abel(x,y): # Li2(x/(1-y))+Li2(y/(1-x))-Li2(xy/((1-x)(1-y)))-Li2(x)-Li2(y)=log(1-x)log(1-y)
    assert N_(x)<1 and N_(y)<1 and N_(x+y)<1
    args=[x/(1-y),y/(1-x),x*y/((1-x)*(1-y)),x,y]
    assert all(N_(t)<1 for t in args)
    return list(zip([1,1,-1,-1,-1],args)), logS(1-x)*logS(1-y)
known={str(canon(b)):Lam, str(canon(-b**4)):M, str(canon(sp.Rational(1,2))):P2/12-lam**2/2, '-1':-P2/12}
def val(x):
    k=str(canon(x)); return known.get(k)
inst=[ ('inversion',inv,(-a**4,),-a**4, -M-8*L**2-P2/6),
 ('Landen',land,(-1/r2,),-1/r2, -Lam-(L-lam/2)**2/2),
 ('duplication',dup,(1/r2,),1/r2, Lam+L**2/2-lam*L/2-lam**2/8+P2/24),
 ('Landen',land,(-a,),-a, -Lam-L**2-P2/24),
 ('inversion',inv,(-b,),-b, Lam+L**2/2-P2/8),
 ('duplication',dup,(b,),b**2, 4*Lam+L**2-P2/4),
 ('reflection',refl,(b**2,),2*b, -4*Lam-3*L**2+2*lam*L+5*P2/12),
 ('Landen',land,(2*b,),-2*a, 4*Lam+L**2-2*lam*L-5*P2/12),
 ('Abel',abel,(-a**2,-2*a),-a**2, -6*Lam-4*L**2+7*P2/24),
 ('reflection',refl,(b,),r2*b, -Lam-L**2+lam*L/2+P2/6),
 ('Landen',land,(r2*b,),-r2, Lam+L**2/2-lam*L/2-P2/6),
 ('Landen',land,(-b/2,),-b/2, -4*Lam-sp.Rational(3,2)*L**2+lam*L-lam**2/2+P2/4),
 ('Abel',abel,(-a,sp.Rational(1,2)),b/(2*r2), -6*Lam-sp.Rational(5,2)*L**2+sp.Rational(3,2)*lam*L-sp.Rational(9,8)*lam**2+11*P2/24),
 ('Landen',land,(-b,),b/r2, -Lam-L**2/2-lam**2/8+P2/8),
 ('Abel',abel,(-a**2,2*b),None, M-10*Lam-L**2-lam*L+7*P2/8),
]
out={}; allok=True
for n,(name,fn,xs,target,printed) in enumerate(inst,1):
    terms,rhs=fn(*xs)
    # numeric residual of the instance itself (110 digits)
    lhs_num=sum(numval(c*Li2(t)) for c,t in terms); res=abs(lhs_num-numval(rhs))
    if target is None:   # line 15: D = Li2(b^2/sqrt2)-Li2(-a^2/sqrt2)
        A=canon(b**2/r2); B=canon(-a**2/r2)
        rest=0; cA=cB=0
        for c,t in terms:
            k=str(canon(t))
            if k==str(A): cA+=c
            elif k==str(B): cB+=c
            else:
                v=val(t); assert v is not None,(n,t); rest+=c*v
        assert cA==1 and cB==-1, (cA,cB)
        got=sp.expand(rhs-rest); key='D15'
    else:
        tk=str(canon(target)); ct=0; rest=0
        for c,t in terms:
            k=str(canon(t))
            if k==tk: ct+=c
            else:
                v=val(t); assert v is not None,(n,t,'needs unknown value'); rest+=c*v
        got=sp.expand((rhs-rest)/ct); known[tk]=got; key=tk
    ok=sp.expand(got-printed)==0
    allok&=ok
    print('(%2d) %-11s x=%s  args=%s'%(n,name,[str(canon(x)) for x in xs],[(str(c),str(canon(t))) for c,t in terms]))
    print('      value = %s   matches printed: %s   instance residual (110 dps) = %s'%(got,ok,mp.nstr(res,3)))
    out[key]=got
# independent check of the obtained values against mpmath Li2 at 110 digits
print('\nvalues vs direct mpmath Li2 (110 dps):')
mp.mp.dps=120
for k,v in known.items():
    if k in ('-1',str(canon(b)),str(canon(-b**4)),str(canon(sp.Rational(1,2)))): 
        if k=='-1' or k=='1/2':
            print('  classical',k,mp.nstr(abs(mp.polylog(2,mp.mpf(sp.N(sp.sympify(k),130)))-numval(v)),3))
        continue
    d=abs(mp.polylog(2,mp.mpf(sp.N(sp.sympify(k),130)))-numval(v.subs({Lam:Li2(b),M:Li2(-b**4)})))
    print('  Li2(%s): %s'%(k,mp.nstr(d,3)))
dv=out['D15']; dd=abs(numval(Li2(b**2/r2)-Li2(-a**2/r2))-numval(dv.subs({Lam:Li2(b),M:Li2(-b**4)})))
print('  D15:',mp.nstr(dd,3))
print('ALL 15 printed values reproduced exactly:',allok)
json.dump({k:sp.srepr(sp.sympify(v)) for k,v in known.items()},open('r3_table_values.json','w'),indent=1)
