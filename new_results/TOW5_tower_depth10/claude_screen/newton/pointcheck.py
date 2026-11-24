# exact-rational point check of Newton identities using OWN rebuilt family (own_F.pkl); no floats, no lambdify
import sys, pickle, random
sys.dont_write_bytecode=True
W="<project>/research/08_theory/TOW5_tower_depth10"
sys.path.insert(0,W)
import sympy as sp
import tow5lib as T5
d,m=T5.d,T5.m; I=sp.I
raw=pickle.load(open("own_F.pkl","rb"))
Fx={k:sp.sympify(v,locals={'d':d,'m':m,'I':I}) for k,v in raw.items()}
def Fv(dp,s,dd,mm):
    return sp.Integer(0) if s>dp//2 else Fx[(dp,s)].subs({d:dd,m:mm},simultaneous=True)
def Lam(delta,s,dd,mm):
    t=0
    for dp in range(delta+1):
        t+=I*(2*dd-1)*(T5.tau1f(delta-dp)*Fv(dp,s,dd+1,mm+1)-T5.tau2f(delta-dp)*Fv(dp,s,dd+1,mm))/(2*dd-dp-1)
    return t
random.seed(5); bad=0; n=0
pts=[(sp.Rational(7),sp.Integer(3)),(sp.Rational(11),sp.Integer(5)),(sp.Rational(19,2),sp.Rational(-3,2)),(sp.Rational(25,3),sp.Integer(1)),(sp.Integer(40),sp.Integer(17))]
for dd,mm in pts:
    for delta in range(1,11):
        for s in range(delta//2+1):
            r=sp.expand(sp.nsimplify(Fv(delta,s,dd,mm)+Fv(delta,s+1,dd,mm)-Lam(delta,s,dd,mm),rational=True))
            n+=1
            if r!=0: bad+=1; print("NEWTON FAIL",delta,s,dd,mm,r)
print("exact points",len(pts),"identities tested",n,"failures",bad)
