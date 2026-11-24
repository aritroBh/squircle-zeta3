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
PA8=sp.sympify(open('P8_A.txt').read()); PC8=sp.sympify(open('P8_CLASS.txt').read())
PB8=sp.sympify(open('P8_B.txt').read()); PZ8=sp.sympify(open('P8_Z.txt').read())
rf=pickle.load(open('refined_gamma_pts.pkl','rb'))
for k,(a1v,(a2r,a4r,a5r)) in rf.items():
    pt={a1s:mp.mpf(a1v),a2s:a2r,a3s:mp.mpf(1),a4s:a4r,a5s:a5r}
    # residuals
    rC4=C4.subs(pt); rQ7=Q7.subs(pt); rP=Q11sp.subs(pt)
    L=a2r*1*a4r*(mp.mpf(a1v)+a2r+1+a4r+a5r)
    A=L*mp.mpf(str(PA8.subs(pt)))/100
    B=L*mp.mpf(str(PB8.subs(pt)))/224
    CL=L*mp.mpf(str(PC8.subs(pt)))/432
    Z=L*mp.mpf(str(PZ8.subs(pt)))/3
    S=mp.mpf(a1v)+a2r+1+a4r+a5r
    print(f"{k} a1={a1v}:")
    print(f"   residuals C4={mp.nstr(abs(rC4),3)} Q7'={mp.nstr(abs(rQ7),3)} P7={mp.nstr(abs(rP),3)}")
    print(f"   A={mp.nstr(A,10)}  B={mp.nstr(B,10)}")
    print(f"   CLASS={mp.nstr(CL,10)}  Z={mp.nstr(Z,10)}")
    print(f"   B/A={mp.nstr(B/A,10)} (75/14={mp.nstr(mp.mpf(75)/14,10)})  Z/A={mp.nstr(Z/A,10)} (-25/3={mp.nstr(-mp.mpf(25)/3,10)})")
    print(f"   A/S^12={mp.nstr(A/S**12,6)} CLASS/S^12={mp.nstr(CL/S**12,6)}")
