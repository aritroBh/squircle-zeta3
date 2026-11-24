import sys, time, pickle, os
from fractions import Fraction as Fr
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import sympy as sp
import tow5lib as T5
d_, m_ = T5.d, T5.m
i = sp.I
t0=time.time()

def qsi_div(x, y):
    def gm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    def gd(a,b):
        den = b[0]**2+b[1]**2; return ((a[0]*b[0]+a[1]*b[1])/den, (a[1]*b[0]-a[0]*b[1])/den)
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2 = (gm(C,C)[0]-2*gm(D,D)[0], gm(C,C)[1]-2*gm(D,D)[1])
    nA = (gm(A,C)[0]-2*gm(B,D)[0], gm(A,C)[1]-2*gm(B,D)[1])
    nB = (gm(B,C)[0]-gm(A,D)[0], gm(B,C)[1]-gm(A,D)[1])
    return (gd(nA, C2), gd(nB, C2))

RHO = pickle.load(open(work+'/data_rho_D10_d7_40.pkl','rb'))
results = {}
for delta in [9, 10]:
    a = pickle.load(open(work+f'/a{delta}_closed.pkl','rb'))
    aa = sp.sympify(a, locals={'d':d_,'m':m_,'I':i})
    bad = 0; n = 0; badex = []
    for dd in RHO:
        if 2*dd-4 < delta: continue
        for key in RHO[dd]:
            p,q = key
            num = RHO[dd][key][delta]; den = RHO[dd][key][0]
            rv = qsi_div(num, den)
            cv = sp.expand(aa.subs({d_: dd, m_: p}))
            exp = (Fr(sp.re(cv)), Fr(sp.im(cv)))
            n += 1
            if rv[1] != (Fr(0),Fr(0)) or rv[0] != exp:
                bad += 1
                if len(badex) < 5: badex.append((dd, key))
    results[delta] = (n, bad, badex)
    print(f"rho({delta})/rho(0) == a_{delta}: {n} cells, {bad} bad  {badex}  ({round(time.time()-t0,1)}s)", flush=True)
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} stepA-rho: a_9 {results[9][0]} cells {results[9][1]} bad; a_10 {results[10][0]} cells {results[10][1]} bad\n")
