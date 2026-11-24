import sympy as sp, mpmath as mp, pickle, os
mp.mp.dps=80
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
def residuals(a1v,a2v,a4v,a5v):
    dd={a1:sp.Float(mp.nstr(a1v,70)),a2:sp.Float(mp.nstr(a2v,70)),a3:1,a4:sp.Float(mp.nstr(a4v,70)),a5:sp.Float(mp.nstr(a5v,70))}
    return max(abs(mp.mpf(str(e.subs(dd)))) for e in (C4,Q7,Q11sp))
def gamma_point(a1v, guess):
    C4s=C4.subs({a1:sp.Rational(a1v),a3:1}); Q7s=Q7.subs({a1:sp.Rational(a1v),a3:1}); Qs=Q11sp.subs({a1:sp.Rational(a1v),a3:1})
    f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    return mp.findroot(f, guess, tol=mp.mpf('1e-60'), maxsteps=1000)
targets=pickle.load(open('gamma_extra.pkl','rb'))
refined=[]
for r in targets:
    a1v=r[0]
    sol=gamma_point(sp.Rational(float(a1v)), (r[1],r[2],r[3]))
    refined.append((a1v,sol[0],sol[1],sol[2]))
for r in refined:
    print(f"  a1={mp.nstr(r[0],6)} resid={mp.nstr(residuals(*r),3)}")
pickle.dump(refined, open('gamma_extra.pkl','wb'))
print("\nproof log tail:")
print('\n'.join(open('proof_div.log').read().split('\n')[-4:]))
