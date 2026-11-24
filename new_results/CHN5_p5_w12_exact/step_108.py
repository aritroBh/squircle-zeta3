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
refpts={}
refpts['p1']=(0.01, refine(0.01,(2.74130851831,2.74173978205,1.04607315237)))
refpts['p2']=(1.7782794100389228, refine(1.7782794100389228,(1.40823387734,2.16140661893,0.362434801222)))
refpts['p3']=(316.22776601683796, refine(316.22776601683796,(197.775002732,306.893601344,113.892966212)))
for k,(a1v,(a2r,a4r,a5r)) in refpts.items():
    print(k, "a1=",a1v, "a2=",mp.nstr(a2r,10),"a4=",mp.nstr(a4r,10),"a5=",mp.nstr(a5r,10))
# evaluate factors
RP=sp.sympify(pickle.load(open('resultant_P7.pkl','rb')))
RQ=sp.sympify(pickle.load(open('resultant_QA.pkl','rb')))
G22=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
flP=sp.factor_list(RP); flQ=sp.factor_list(RQ)
facsP=[(sp.expand(f),m) for f,m in flP[1]]
facsQ=[(sp.expand(f),m) for f,m in flQ[1]]
def ev(P, a1v, s):
    return P.subs({a1s:mp.mpf(a1v),a2s:s[0],a4s:s[1]})
print("\nG22 and factor values (should be ~0 on the containing factor):")
for k,(a1v,s) in refpts.items():
    print(f"  {k}: G22={mp.nstr(abs(mp.mpf(str(ev(G22,a1v,s)))),3)}")
    for i,(F,m) in enumerate(facsP):
        print(f"      ResP fac{i} (deg{sp.Poly(F).total_degree()}): {mp.nstr(abs(mp.mpf(str(ev(F,a1v,s)))),3)}")
    for i,(F,m) in enumerate(facsQ):
        print(f"      ResQ fac{i} (deg{sp.Poly(F).total_degree()}): {mp.nstr(abs(mp.mpf(str(ev(F,a1v,s)))),3)}")
