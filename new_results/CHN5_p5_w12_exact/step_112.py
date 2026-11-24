import os, sys, pickle, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp
mp.mp.dps=100
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
RQ=sp.sympify(pickle.load(open('resultant_QA.pkl','rb')))
G30=sp.expand(sp.factor_list(RQ)[1][2][0])
G30l=sp.lambdify((a1s,a2s,a4s), G30, 'mpmath')
res=pickle.load(open('gamma_trace2.pkl','rb'))
pos=0; neg=0; near0=mp.mpf('1e100'); cnt=0
for a1v in sorted(res.keys()):
    for s in res[a1v]:
        if s['a5']>mp.mpf('1e-25'):
            v=G30l(mp.mpf(a1v),s['a2'],s['a4'])
            cnt+=1
            if v>0: pos+=1
            else: neg+=1
            if abs(v)<near0: near0=abs(v)
print(f"G30 (QA-resultant condition) at {cnt} traced Gamma pts: positive={pos}, negative={neg}, min|G30|={mp.nstr(near0,3)}")
# also at the refined high-precision points
rf=pickle.load(open('refined_gamma_pts.pkl','rb')) if os.path.exists('refined_gamma_pts.pkl') else None
