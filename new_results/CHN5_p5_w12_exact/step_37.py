import os, sys, time
import numpy as np, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp
import survey_lib as SL
mp.mp.dps=60
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1s,a2s,a3s,a4s,a5s][i]**e[i] for i in range(5)) for e,c in SL.Q11.items())
def gamma_at(a1v, guesses):
    """find all Gamma points at given a1 (a3=1): solve C4=Q7'=Q11=0 in (a2,a4,a5)."""
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    sols=[]
    for g in guesses:
        try:
            r=mp.findroot(f, g, tol=mp.mpf('1e-45'), maxsteps=100)
            a2r,a4r,a5r=r[0],r[1],r[2]
            if a2r>0 and a4r>0 and a5r>0:
                key=tuple(mp.nstr(x,30) for x in (a2r,a4r,a5r))
                if not any(key==k for k,_ in sols):
                    sols.append((key,(a2r,a4r,a5r)))
        except Exception:
            pass
    return [s for _,s in sols]
# generate guesses at a1 from V10 scan over a2 grid
def guesses_from_scan(a1v, n=120, lo=-4, hi=4):
    gs=[]
    for a2v in np.logspace(lo,hi,n):
        for (a4v,a5v) in SL.v10_points(a1v, a2v):
            gs.append((mp.mpf(a2v),mp.mpf(a4v),mp.mpf(a5v)))
    return gs
t0=time.time()
for a1v in [0.01, 1.0044068, 2.0]:
    gs=guesses_from_scan(a1v, n=60, lo=-2, hi=1.5)
    sols=gamma_at(a1v, gs)
    print(f"a1={a1v}: {len(gs)} guesses -> {len(sols)} Gamma pts")
    for (a2r,a4r,a5r) in sols:
        print(f"   a2={mp.nstr(a2r,12)} a4={mp.nstr(a4r,12)} a5={mp.nstr(a5r,12)}")
print(f"[{time.time()-t0:.0f}s]")
