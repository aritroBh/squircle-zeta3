import sympy as sp, mpmath as mp, pickle
mp.mp.dps=70
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
def gamma_point(a1v, guess, a2free=True):
    """solve C4=Q7'=Q11=0 at a1=a1v. unknowns (a2,a4,a5)."""
    C4s=C4.subs({a1:sp.Rational(a1v),a3:1}); Q7s=Q7.subs({a1:sp.Rational(a1v),a3:1}); Qs=Q11sp.subs({a1:sp.Rational(a1v),a3:1})
    f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    sol=mp.findroot(f, guess, tol=mp.mpf('1e-55'), maxsteps=200)
    return sol
sol=gamma_point(sp.Rational(1,100), (mp.mpf('2.74'),mp.mpf('2.74'),mp.mpf('1.046')))
a2r,a4r,a5r=sol[0],sol[1],sol[2]
print("Gamma point at a1=1/100:")
print(f"  a2={mp.nstr(a2r,40)}")
print(f"  a4={mp.nstr(a4r,40)}")
print(f"  a5={mp.nstr(a5r,40)}")
# verify residuals
pt=[mp.mpf('0.01'),a2r,mp.mpf(1),a4r,a5r]
for nm,e in [('C4',C4),('Q7',Q7),('Q11',Q11sp)]:
    v=e.subs({a1:pt[0],a2:sp.Float(mp.nstr(a2r,60)),a3:1,a4:sp.Float(mp.nstr(a4r,60)),a5:sp.Float(mp.nstr(a5r,60))})
    print(f"  {nm} residual = {mp.nstr(abs(mp.mpf(str(v))),5)}")
pickle.dump({'a1':mp.mpf('0.01'),'a2':a2r,'a4':a4r,'a5':a5r}, open('gamma_pt1.pkl','wb'))
print("saved gamma_pt1.pkl")
