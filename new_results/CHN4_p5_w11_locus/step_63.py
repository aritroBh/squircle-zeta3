import sympy as sp, mpmath as mp, pickle
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
def residuals(a1v,a2v,a4v,a5v):
    d={a1:sp.Float(mp.nstr(a1v,50)),a2:sp.Float(mp.nstr(a2v,50)),a3:1,a4:sp.Float(mp.nstr(a4v,50)),a5:sp.Float(mp.nstr(a5v,50))}
    return [abs(mp.mpf(str(e.subs(d)))) for e in (C4,Q7,Q11sp)]
def gamma_point_multi(a1v, guesses):
    C4s=C4.subs({a1:sp.Rational(a1v),a3:1}); Q7s=Q7.subs({a1:sp.Rational(a1v),a3:1}); Qs=Q11sp.subs({a1:sp.Rational(a1v),a3:1})
    f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    sols=[]
    for g in guesses:
        try:
            s=mp.findroot(f, g, tol=mp.mpf('1e-45'), maxsteps=300)
            pt=(s[0],s[1],s[2])
            if all(x>0 for x in pt):
                rr=residuals(mp.mpf(a1v),*pt)
                if max(rr)<mp.mpf('1e-30') and not any(all(abs(pt[i]-q[i])<mp.mpf('1e-20') for i in range(3)) for q in sols):
                    sols.append(pt)
        except Exception: pass
    return sols
# verify up branch residuals
up=pickle.load(open('gamma_up.pkl','rb'))
print("up-branch residuals (max of C4,Q7,Q11):")
for r in up:
    rr=residuals(*r)
    print(f"  a1={mp.nstr(r[0],4)}: max resid={mp.nstr(max(rr),4)}  a2={mp.nstr(r[1],6)} a4={mp.nstr(r[2],6)} a5={mp.nstr(r[3],6)}")
# independently probe several a1 with a grid of guesses
print("\nindependent probes:")
import numpy as np
for a1v in [0.5, 1.0, 2.0, 5.0]:
    guesses=[(mp.mpf(a2),mp.mpf(a4),mp.mpf(a5)) for a2 in np.linspace(0.3,6,8) for a4 in np.linspace(0.3,6,8) for a5 in np.linspace(0.2,3,5)]
    sols=gamma_point_multi(sp.Rational(a1v), guesses)
    print(f"  a1={a1v}: {len(sols)} Gamma pts: "+"; ".join(f"(a2={mp.nstr(s[0],5)},a4={mp.nstr(s[1],5)},a5={mp.nstr(s[2],5)})" for s in sols[:4]))
