import sys, time, pickle
sys.dont_write_bytecode=True
W="<project>/research/08_theory/TOW5_tower_depth10"
sys.path.insert(0,W)
import sympy as sp
import tow5lib as T5
d,m=T5.d,T5.m
I=sp.I
t1=T5.tau1f; t2=T5.tau2f
# own rebuild: F(delta,0)=arho(delta) (all delta<=10), recursion own implementation, no pickled F
F={}
def Fv(dp,s):
    return F[(dp,s)] if s<=dp//2 else sp.Integer(0)
def Lam(delta,s):
    tot=0
    for dp in range(delta+1):
        fA=Fv(dp,s).subs({d:d+1,m:m+1},simultaneous=True)
        fB=Fv(dp,s).subs({d:d+1,m:m},simultaneous=True)
        tot+= I*(2*d-1)*(t1(delta-dp)*fA-t2(delta-dp)*fB)/(2*d-dp-1)
    return tot
t0=time.time()
for delta in range(0,11):
    F[(delta,0)]=sp.cancel(T5.arho(delta))
    for s in range(0,delta//2+1):
        nxt=sp.cancel(sp.together(Lam(delta,s)-Fv(delta,s)))
        if s+1<=delta//2: F[(delta,s+1)]=nxt
        else:
            print(f"closing F({delta},{s+1}) == 0 (cancel'd rational function is 0):",nxt==0, f"{time.time()-t0:.0f}s",flush=True)
# compare with pickled forced family
raw=pickle.load(open(W+"/forced_F910.pkl","rb"))
bad=0;n=0
for k,v in raw.items():
    if isinstance(k,tuple) and k[0]!='closed' and k in F:
        n+=1
        diff=sp.cancel(sp.sympify(v,locals={'d':d,'m':m,'I':I})-F[k])
        if diff!=0: bad+=1; print("MISMATCH",k)
print("compared",n,"pickled cells vs own rebuild; mismatches",bad)
print("closed flags in pickle:",{k:v for k,v in raw.items() if k[0]=='closed'} )
pickle.dump({k:str(v) for k,v in F.items()},open("own_F.pkl","wb"))
# sample-point closing check, exact rationals, F(9,5), F(10,6) evaluated by recursion at points via the cancelled forms
