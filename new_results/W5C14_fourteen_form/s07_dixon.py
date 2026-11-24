"""Exact solve y . R_{S,C} = Phi_s(Delta)_C per sector by Dixon lifting, REUSING W5C2's trimmed systems
(S = independent rows, C = pivot columns; read-only from W5C2_exact_certificate/work).  Usage: s07_dixon.py pp|pm|mp|mm"""
import sys, os, time, pickle
sys.dont_write_bytecode=True
TH='<project>/research/08_theory'
W5C2=TH+'/W5C2_exact_certificate'
sys.path.insert(0, W5C2)
import numpy as np
from fractions import Fraction as F
from math import lcm
import w5c2_dixon as dx
P1=16777213; P2=16777199
def log(s): print(time.strftime('%H:%M:%S ')+s, flush=True)
secname=sys.argv[1]; sec={'pp':(1,1),'pm':(1,-1),'mp':(-1,1),'mm':(-1,-1)}[secname]
t0=time.monotonic()
d2=pickle.load(open(f'{W5C2}/work/step02_{secname}.pkl','rb'))
S,C=d2['S'],d2['C']
rows=pickle.load(open(f'{W5C2}/work/feed_{secname}.pkl','rb'))
r=len(S); Cset=sorted(C); cpos={c:k for k,c in enumerate(Cset)}
A=np.zeros((r,r),np.int64)
for j,si in enumerate(S):
    for c,a in rows[si][2]:
        k=cpos.get(c)
        if k is not None: A[k,j]=a
log(f'{secname}: A {r}x{r} filled maxentry={int(np.abs(A).max())} ({time.monotonic()-t0:.0f}s)')
del rows
folds=pickle.load(open('delta_folds.pkl','rb'))
Bs=np.zeros((r,2),np.int64); Ls={}
for jj,comp in enumerate(('rat','sq')):
    D=folds[(comp,sec)]
    DC={c:v for c,v in D.items() if c in cpos}
    L=1
    for v in DC.values(): L=lcm(L,v.denominator)
    Ls[comp]=L
    for c,v in DC.items():
        x=int(v*L); assert abs(x)<2**62; Bs[cpos[c],jj]=x
log(f'{secname}: B built Ls={Ls} maxB={int(np.abs(Bs).max())}')
Y,steps,M=dx.dixon_solve(A,Bs,P1,ckpt_path=f'dixon_ckpt_{secname}.pkl')
assert Y is not None
log(f'{secname}: CONVERGED steps={steps} M_digits={len(str(M))} ({time.monotonic()-t0:.0f}s)')
ok=dx.check_mod(Y,A,Bs,P2); log(f'{secname}: cross-check A.Y==B mod p2: {ok}'); assert ok
out={'sec':secname,'steps':steps,'Ls':Ls}
for jj,comp in enumerate(('rat','sq')):
    y=[F(Y[i][jj][0],Y[i][jj][1])/Ls[comp] for i in range(r)]
    out['y_'+comp]=y
    nz=sum(1 for a in y if a)
    hn=max((len(str(abs(a.numerator))) for a in y),default=1); hd=max((len(str(a.denominator)) for a in y),default=1)
    log(f'{secname}: y_{comp} nonzero={nz}/{r} max_num_digits={hn} max_den_digits={hd}')
pickle.dump(out,open(f'yexact_{secname}.pkl','wb'),protocol=4)
log(f'{secname}: saved yexact_{secname}.pkl total {time.monotonic()-t0:.0f}s')
