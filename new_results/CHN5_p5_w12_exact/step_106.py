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
    r=mp.findroot(f, tuple(mp.mpf(x) for x in guess), tol=mp.mpf('1e-120'), maxsteps=400)
    return r[0],r[1],r[2]
# refine 3 points
pts={}
pts['p1']=refine(0.01, (2.7413085,2.7417398,1.0460732))
pts['p2']=refine(1.7782794100389228, (3.9535052,2.5565539,3.6050562))  # one of the a1=1.78 points
pts['p3']=refine(316.22776601683796, None if False else (9.1269,9.1540,4.846))
for k,(a2r,a4r,a5r) in pts.items():
    print(k, mp.nstr(a2r,12), mp.nstr(a4r,12), mp.nstr(a5r,12))
pickle.dump({'p1':pts['p1'],'p2':pts['p2'],'p3':pts['p3'],
             'a1':{'p1':0.01,'p2':1.7782794100389228,'p3':316.22776601683796}}, open('refined_gamma_pts.pkl','wb'))
