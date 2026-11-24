import os, pickle, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
mp.mp.dps=50
import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
RP=sp.sympify(pickle.load(open('resultant_P7.pkl','rb')))
RQ=sp.sympify(pickle.load(open('resultant_QA.pkl','rb')))
flP=sp.factor_list(RP); flQ=sp.factor_list(RQ)
facsP=[(sp.expand(f),m) for f,m in flP[1]]
facsQ=[(sp.expand(f),m) for f,m in flQ[1]]
G22=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
res=pickle.load(open('gamma_trace2.pkl','rb'))
testpts=[]
for a1v in [0.01, 1.7782794100389228, 316.22776601683796]:
    for s in res[a1v]:
        if s['a5']>mp.mpf('1e-25'):
            testpts.append((s['a1'],s['a2'],s['a4']))
def ev(P, pt):
    return mp.mpf(str(P.subs({a1:mp.mpf(pt[0]),a2:mp.mpf(pt[1]),a4:mp.mpf(pt[2])})))
print(f"testing at {len(testpts)} traced Gamma points")
for i,(F,m) in enumerate(facsP):
    vals=[mp.nstr(abs(ev(F,p)),3) for p in testpts]
    print(f"  ResP factor {i} (deg {sp.Poly(F).total_degree()}, mult {m}): {vals}")
for i,(F,m) in enumerate(facsQ):
    vals=[mp.nstr(abs(ev(F,p)),3) for p in testpts]
    print(f"  ResQ factor {i} (deg {sp.Poly(F).total_degree()}, mult {m}): {vals}")
vals=[mp.nstr(abs(ev(G22,p)),3) for p in testpts]
print(f"  G22: {vals}")
