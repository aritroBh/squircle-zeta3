import os, sys, time
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp, mpmath as mp
import survey_lib as SL
mp.mp.dps=80
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1s,a2s,a3s,a4s,a5s][i]**e[i] for i in range(5)) for e,c in SL.Q11.items())
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
# mpf-based V10 solver at given (a1,a2): roots of G in a4 via mp.polyroots, then common a5
def v10_mpf(a1v,a2v):
    Gp=sp.Poly(G.subs({a1s:mp.mpf(a1v),a2s:mp.mpf(a2v)}), a4s)
    coeffs=[mp.mpf(str(c)) for c in Gp.all_coeffs()]
    try: roots=mp.polyroots(coeffs, maxsteps=200, extraprec=50)
    except Exception as e: return ('ERR',str(e))
    out=[]
    for r in roots:
        if abs(mp.im(r))<mp.mpf('1e-30') and mp.re(r)>0:
            a4v=mp.re(r)
            # common a5 of C4,Q7' at (a1,a2,a4)
            C4p=sp.Poly(C4.subs({a1s:mp.mpf(a1v),a2s:mp.mpf(a2v),a4s:a4v,a3s:1}), a5s)
            Q7p=sp.Poly(Q7.subs({a1s:mp.mpf(a1v),a2s:mp.mpf(a2v),a4s:a4v,a3s:1}), a5s)
            cc=[mp.mpf(str(c)) for c in C4p.all_coeffs()]
            try: a5roots=mp.polyroots(cc, maxsteps=200, extraprec=50)
            except Exception: continue
            for a5r in a5roots:
                if abs(mp.im(a5r))<mp.mpf('1e-25') and mp.re(a5r)>0:
                    qv=[mp.mpf(str(c)) for c in Q7p.all_coeffs()]
                    val=mp.polyval(qv, a5r)
                    if abs(val)<mp.mpf('1e-20'):
                        out.append((a4v,mp.re(a5r)))
    return out
# check at a1=38311 (where trace found 0) and a1=14678 (found 1)
for a1v in [14678, 38311, 100000]:
    t0=time.time()
    # scan a2 for V10 points and Q11 sign
    import numpy as np
    found=0; npos=0
    q11vals=[]
    for a2v in np.logspace(-3,5,60):
        r=v10_mpf(a1v,a2v)
        if r and r[0]!='ERR':
            for (a4v,a5v) in r:
                npos+=1
                q=SL.q11_np([a1v,a2v,1,float(a4v),float(a5v)])
                q11vals.append(q)
    sgns=set('+' if q>0 else '-' for q in q11vals)
    print(f"a1={a1v}: {npos} positive V10 sheet-pts scanned, Q11 signs present: {sgns} [{time.time()-t0:.0f}s]")
