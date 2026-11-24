import sympy as sp, mpmath as mp, pickle
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
def residuals(a1v,a2v,a4v,a5v):
    d={a1:sp.Float(mp.nstr(a1v,50)),a2:sp.Float(mp.nstr(a2v,50)),a3:1,a4:sp.Float(mp.nstr(a4v,50)),a5:sp.Float(mp.nstr(a5v,50))}
    return [abs(mp.mpf(str(e.subs(d)))) for e in (C4,Q7,Q11sp)]
up=pickle.load(open('gamma_up.pkl','rb'))
print("up-branch residuals:")
for r in up:
    print(f"  a1={mp.nstr(r[0],4)}: max resid={mp.nstr(max(residuals(*r)),4)}")
# trace DOWN from a1=0.01 (g0)
def gamma_point(a1v, guess):
    C4s=C4.subs({a1:sp.Rational(a1v),a3:1}); Q7s=Q7.subs({a1:sp.Rational(a1v),a3:1}); Qs=Q11sp.subs({a1:sp.Rational(a1v),a3:1})
    f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    return mp.findroot(f, guess, tol=mp.mpf('1e-45'), maxsteps=300)
g0=(mp.mpf('2.7413'),mp.mpf('2.7417'),mp.mpf('1.0461'))
down=[]; cur=g0
for a1v in [0.005,0.002,0.001,0.0005]:
    try:
        sol=gamma_point(sp.Rational(a1v), cur); cur=(sol[0],sol[1],sol[2])
        if all(x>0 for x in cur): down.append((mp.mpf(a1v),cur[0],cur[1],cur[2]))
    except Exception as ex:
        print(f"down failed at a1={a1v}"); break
print("\ndown branch:")
for r in down:
    print(f"  a1={mp.nstr(r[0],6)}: resid={mp.nstr(max(residuals(*r)),4)} a2={mp.nstr(r[1],6)} a4={mp.nstr(r[2],6)} a5={mp.nstr(r[3],6)}")
pickle.dump(down, open('gamma_down.pkl','wb'))
