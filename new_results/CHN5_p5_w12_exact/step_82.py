import os
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W)
trace2 = r'''
import sys, time, pickle
import sympy as sp, mpmath as mp
import numpy as np
sys.path.insert(0,'.')
import survey_lib as SL
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1s,a2s,a3s,a4s,a5s][i]**e[i] for i in range(5)) for e,c in SL.Q11.items())
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
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
def q11_mpf(pt): return mpoly(SL.Q11, pt)
def v10_mpf(a1v,a2v):
    Gp=sp.Poly(G.subs({a1s:mp.mpf(a1v),a2s:mp.mpf(a2v)}), a4s)
    try: roots=mp.polyroots([mp.mpf(str(c)) for c in Gp.all_coeffs()], maxsteps=300, extraprec=60)
    except Exception: return []
    out=[]
    for r in roots:
        if abs(mp.im(r))<mp.mpf('1e-15')*max(1,abs(mp.re(r))) and mp.re(r)>0:
            a4v=mp.re(r)
            C4p=sp.Poly(C4.subs({a1s:mp.mpf(a1v),a2s:mp.mpf(a2v),a4s:a4v,a3s:1}), a5s)
            Q7p=sp.Poly(Q7.subs({a1s:mp.mpf(a1v),a2s:mp.mpf(a2v),a4s:a4v,a3s:1}), a5s)
            try: a5roots=mp.polyroots([mp.mpf(str(c)) for c in C4p.all_coeffs()], maxsteps=300, extraprec=60)
            except Exception: continue
            for a5r in a5roots:
                if abs(mp.im(a5r))<mp.mpf('1e-12')*max(1,abs(mp.re(a5r))) and mp.re(a5r)>0:
                    if abs(mp.polyval([mp.mpf(str(c)) for c in Q7p.all_coeffs()], a5r))<mp.mpf('1e-10'):
                        out.append((a4v,mp.re(a5r)))
    return out
def gamma_at_mpf(a1v, cont, n=80, lo=-4, hi=6):
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    guesses=list(cont)
    prev=None
    for a2v in np.logspace(lo,hi,n):
        cur=[]
        for (a4v,a5v) in v10_mpf(a1v,a2v):
            cur.append((a4v,a5v,q11_mpf([a1v,a2v,1,a4v,a5v])))
        if prev:
            a2p,prevl=prev
            for (a4v,a5v,q) in cur:
                if not prevl: continue
                j=min(range(len(prevl)), key=lambda k:(prevl[k][0]-a4v)**2+(prevl[k][1]-a5v)**2)
                a4p,a5p,qp=prevl[j]
                if qp*q<0: guesses.append(((a2p+a2v)/2,(a4p+a4v)/2,(a5p+a5v)/2))
        prev=(a2v,cur)
    sols=[]
    for g in guesses:
        try:
            r=mp.findroot(f, tuple(mp.mpf(x) for x in g), tol=mp.mpf('1e-40'), maxsteps=150)
            a2r,a4r,a5r=r[0],r[1],r[2]
            if a2r>0 and a4r>0 and a5r>0:
                key=tuple(mp.nstr(x,18) for x in (a2r,a4r,a5r))
                if not any(key==k for k,_ in sols): sols.append((key,(a2r,a4r,a5r)))
        except Exception: pass
    out=[]
    for _,(a2r,a4r,a5r) in sols:
        pt=[mp.mpf(a1v),a2r,mp.mpf(1),a4r,a5r]
        qa=mpoly(QA,pt); qc=mpoly(QC,pt)
        out.append({'a1':mp.mpf(a1v),'a2':a2r,'a4':a4r,'a5':a5r,'QA':qa,'QC':qc})
    return out
schedule=list(np.logspace(-7,6,53))
results={}
lg=open('gamma_trace2.log','w')
prev=[]; t0=time.time()
for a1v in schedule:
    sols=gamma_at_mpf(a1v, prev)
    results[a1v]=sols
    prev=[(s['a2'],s['a4'],s['a5']) for s in sols]
    sgA=sorted(set('+' if s['QA']>0 else '-' for s in sols))
    sgC=sorted(set('+' if s['QC']>0 else '-' for s in sols))
    lg.write(f"a1={a1v:.6g}: {len(sols)} pts signQA={sgA} signQC={sgC} [{time.time()-t0:.0f}s]\n"); lg.flush()
    pickle.dump(results, open('gamma_trace2_partial.pkl','wb'))
pickle.dump(results, open('gamma_trace2.pkl','wb'))
lg.write("DONE\n"); lg.close()
print("trace2 done")
'''
open(W+'/run_gamma_trace2.py','w').write(trace2)
print("run_gamma_trace2.py written")
import subprocess, time
proc=subprocess.Popen(['nice','-n','10','<project>/.venv/bin/python','run_gamma_trace2.py'],
                      cwd=W, stdout=open(W+'/run_gamma_trace2.out','w'), stderr=subprocess.STDOUT)
open(W+'/run_gamma_trace2.pid','w').write(str(proc.pid))
print("launched trace2 pid", proc.pid, "at", time.strftime('%H:%M:%S'))
