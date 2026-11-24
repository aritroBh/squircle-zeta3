import os
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
trace = r'''
import sys, time, pickle
import numpy as np, mpmath as mp
sys.path.insert(0,'.')
import sympy as sp
import survey_lib as SL
from fractions import Fraction as F
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1s,a2s,a3s,a4s,a5s][i]**e[i] for i in range(5)) for e,c in SL.Q11.items())
Q=pickle.load(open('Q8_quotients.pkl','rb'))
mp.mp.dps=60
def mpoly(Qd, pt):
    tot=mp.mpf(0)
    for e,c in Qd.items():
        pw=mp.mpf(1)
        for x,ei in zip(pt,e):
            if ei: pw*=x**ei
        tot+=mp.mpf(c.numerator)/c.denominator*pw
    return tot
QA=Q['A']; QC=Q['CLASS']
def gamma_points(a1v, cont_guesses):
    mp.mp.dps=35
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    guesses=list(cont_guesses)
    for a2v in np.logspace(-3,3,80):
        for (a4v,a5v) in SL.v10_points(a1v, a2v):
            guesses.append((mp.mpf(a2v),mp.mpf(a4v),mp.mpf(a5v)))
    cand=[]
    for g in guesses:
        try:
            r=mp.findroot(f, g, tol=mp.mpf('1e-28'), maxsteps=80)
            a2r,a4r,a5r=r[0],r[1],r[2]
            if a2r>mp.mpf('1e-9') and a4r>mp.mpf('1e-9') and a5r>mp.mpf('1e-9'):
                key=tuple(mp.nstr(x,18) for x in (a2r,a4r,a5r))
                if not any(key==k for k,_ in cand): cand.append((key,(a2r,a4r,a5r)))
        except Exception: pass
    # refine at dps=60
    mp.mp.dps=60
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    out=[]
    for _,(a2r,a4r,a5r) in cand:
        try:
            r=mp.findroot(f, (a2r,a4r,a5r), tol=mp.mpf('1e-50'), maxsteps=200)
            a2r,a4r,a5r=r[0],r[1],r[2]
            pt=[mp.mpf(a1v),a2r,mp.mpf(1),a4r,a5r]
            res=max(abs(fi(a2r,a4r,a5r)) for fi in f)
            qa=mpoly(QA,pt); qc=mpoly(QC,pt)
            out.append({'a1':mp.mpf(a1v),'a2':a2r,'a4':a4r,'a5':a5r,'QA':qa,'QC':qc,'res':res})
        except Exception: pass
    return out
schedule=list(np.logspace(-5,5,121))
results={}
lg=open('gamma_trace.log','w')
prev=[]
t0=time.time()
for a1v in schedule:
    sols=gamma_points(a1v, prev)
    results[a1v]=sols
    prev=[(s['a2'],s['a4'],s['a5']) for s in sols]
    sgA=set('+' if s['QA']>0 else '-' for s in sols)
    sgC=set('+' if s['QC']>0 else '-' for s in sols)
    lg.write(f"a1={a1v:.6g}: {len(sols)} Gamma pts; signQA={sorted(sgA)} signQC={sorted(sgC)} [{time.time()-t0:.0f}s]\n"); lg.flush()
    pickle.dump(results, open('gamma_trace_partial.pkl','wb'))
pickle.dump(results, open('gamma_trace.pkl','wb'))
lg.write("DONE\n"); lg.close()
print("trace done")
'''
open(W+'/run_gamma_trace.py','w').write(trace)
print("run_gamma_trace.py written")
