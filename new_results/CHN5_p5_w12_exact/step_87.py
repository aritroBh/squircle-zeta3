import os, sys
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp, mpmath as mp
import survey_lib as SL
import pickle
mp.mp.dps=60
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1s,a2s,a3s,a4s,a5s][i]**e[i] for i in range(5)) for e,c in SL.Q11.items())
Q=pickle.load(open('Q8_quotients.pkl','rb'))
def mpoly(Qd, pt):
    tot=mp.mpf(0)
    for e,c in Qd.items():
        pw=mp.mpf(1)
        for x,ei in zip(pt,e):
            if ei: pw*=x**ei
        tot+=mp.mpf(c.numerator)/c.denominator*pw
    return tot
QA=Q['A']; QC=Q['CLASS']
def solve_near(a1v, guess):
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    try:
        r=mp.findroot(f, tuple(mp.mpf(x) for x in guess), tol=mp.mpf('1e-45'), maxsteps=300)
        return r[0],r[1],r[2]
    except Exception as e:
        return None
# probe near-stratum interior at a1=1+-0.02, near a2=a4=8241, 249, 1.758, small a5
for a1v in [0.98, 1.02]:
    print(f"--- a1={a1v} ---")
    for a2g in [1.758, 249.46, 8241.15]:
        for a5g in [0.001, 0.01]:
            s=solve_near(a1v,(a2g,a2g,a5g))
            if s:
                a2r,a4r,a5r=s
                if a2r>0 and a4r>0 and a5r>1e-8:
                    pt=[mp.mpf(a1v),a2r,mp.mpf(1),a4r,a5r]
                    qa=mpoly(QA,pt); qc=mpoly(QC,pt)
                    print(f"  guess a2={a2g}: a2={mp.nstr(a2r,10)} a4={mp.nstr(a4r,10)} a5={mp.nstr(a5r,6)} QA={'+' if qa>0 else '-'} QC={'+' if qc>0 else '-'} |QA|={mp.nstr(abs(qa),3)}")
