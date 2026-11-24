import os, sys, pickle, mpmath as mp
base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
os.chdir(W); sys.path.insert(0,W)
import sympy as sp
import survey_lib as SL
mp.mp.dps=150
a1s,a2s,a3s,a4s,a5s=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1s,a2s,a3s,a4s,a5s][i]**e[i] for i in range(5)) for e,c in SL.Q11.items())
def refine(a1v, guess):
    C4s=C4.subs({a1s:mp.mpf(a1v),a3s:1}); Q7s=Q7.subs({a1s:mp.mpf(a1v),a3s:1}); Qs=Q11sp.subs({a1s:mp.mpf(a1v),a3s:1})
    f=[sp.lambdify((a2s,a4s,a5s), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    r=mp.findroot(f, tuple(mp.mpf(x) for x in guess), tol=mp.mpf('1e-110'), maxsteps=400)
    return r[0],r[1],r[2]
rf={'p1':(0.01, refine(0.01,(2.74130851831,2.74173978205,1.04607315237))),
    'p2':(1.7782794100389228, refine(1.7782794100389228,(1.40823387734,2.16140661893,0.362434801222))),
    'p3':(316.22776601683796, refine(316.22776601683796,(197.775002732,306.893601344,113.892966212)))}
pickle.dump(rf, open('refined_gamma_pts.pkl','wb'))
PA8=sp.sympify(open('P8_A.txt').read()); PC8=sp.sympify(open('P8_CLASS.txt').read())
PB8=sp.sympify(open('P8_B.txt').read()); PZ8=sp.sympify(open('P8_Z.txt').read())
for k,(a1v,(a2r,a4r,a5r)) in rf.items():
    pt={a1s:mp.mpf(a1v),a2s:a2r,a3s:mp.mpf(1),a4s:a4r,a5s:a5r}
    rC4=C4.subs(pt); rQ7=Q7.subs(pt); rP=Q11sp.subs(pt)
    L=a2r*a4r*(mp.mpf(a1v)+a2r+1+a4r+a5r)
    A=L*mp.mpf(str(PA8.subs(pt)))/100; B=L*mp.mpf(str(PB8.subs(pt)))/224
    CL=L*mp.mpf(str(PC8.subs(pt)))/432; Z=L*mp.mpf(str(PZ8.subs(pt)))/3
    S=mp.mpf(a1v)+a2r+1+a4r+a5r
    print(f"{k} a1={a1v}: residuals(C4,Q7',P7)=({mp.nstr(abs(rC4),2)},{mp.nstr(abs(rQ7),2)},{mp.nstr(abs(rP),2)})")
    print(f"   A={mp.nstr(A,12)}  sign={'+' if A>0 else '-'}   A/S^12={mp.nstr(A/S**12,6)}")
    print(f"   CLASS={mp.nstr(CL,12)} sign={'+' if CL>0 else '-'}  B/A={mp.nstr(B/A,8)} Z/A={mp.nstr(Z/A,8)}")
