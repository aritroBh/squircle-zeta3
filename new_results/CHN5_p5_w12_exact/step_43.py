import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp
import survey_lib as SL
mp.mp.dps=50
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1s,a2s,a3s,a4s,a5s][i]**e[i] for i in range(5)) for e,c in SL.Q11.items())

def scan_sheets(a1v, n=160, lo=-4, hi=4):
    """return list of (a2lo,a2hi, guess(a2,a4,a5)) brackets where Q11 changes sign on a V10 sheet."""
    a2grid=np.logspace(lo,hi,n)
    prev=None  # list of (a4,a5,q11) at prev a2
    brackets=[]
    for i,a2v in enumerate(a2grid):
        cur=[]
        for (a4v,a5v) in SL.v10_points(a1v, a2v):
            q=SL.q11_np([a1v,a2v,1,a4v,a5v])
            cur.append((a4v,a5v,q))
        if prev is not None:
            a2p, prevl = prev
            for (a4v,a5v,q) in cur:
                # nearest prev point in (a4,a5)
                if not prevl: continue
                j=min(range(len(prevl)), key=lambda k:(prevl[k][0]-a4v)**2+(prevl[k][1]-a5v)**2)
                a4p,a5p,qp=prevl[j]
                if qp*q<0 and abs(a4p-a4v)<0.5*max(1,abs(a4v)):
                    brackets.append((a2p,a2v,(a4p+a4v)/2,(a5p+a5v)/2))
        prev=(a2v,cur)
    return brackets

def gamma_at_fast(a1v, dps=50):
    mp.mp.dps=dps
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    sols=[]
    for (a2lo,a2hi,a4g,a5g) in scan_sheets(a1v):
        try:
            r=mp.findroot(f, (mp.mpf((a2lo+a2hi)/2),mp.mpf(a4g),mp.mpf(a5g)), tol=mp.mpf('1e-40'), maxsteps=100)
            a2r,a4r,a5r=r[0],r[1],r[2]
            if a2r>1e-12 and a4r>1e-12 and a5r>1e-12:
                key=tuple(mp.nstr(x,25) for x in (a2r,a4r,a5r))
                if not any(key==k for k,_ in sols): sols.append((key,(a2r,a4r,a5r)))
        except Exception: pass
    return [s for _,s in sols]

t0=time.time()
for a1v in [0.01, 1.0044068, 2.0]:
    sols=gamma_at_fast(a1v)
    print(f"a1={a1v}: {len(sols)} Gamma pts")
    for (a2r,a4r,a5r) in sols:
        print(f"   a2={mp.nstr(a2r,12)} a4={mp.nstr(a4r,12)} a5={mp.nstr(a5r,12)}")
print(f"[{time.time()-t0:.0f}s]")
