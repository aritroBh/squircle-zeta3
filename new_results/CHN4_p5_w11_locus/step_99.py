import sympy as sp, mpmath as mp, pickle
mp.mp.dps=50
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
# boundary a1=0 system in (a2,a4,a5), a3=1
C4b=C4.subs({a1:0,a3:1}); Q7b=Q7.subs({a1:0,a3:1}); Qb=Q11sp.subs({a1:0,a3:1})
f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4b,Q7b,Qb)]
sol=mp.findroot(f, (mp.mpf('2.667'),mp.mpf('2.667'),mp.mpf('1.0')), tol=mp.mpf('1e-40'), maxsteps=400)
print("a1=0 boundary limit point:")
print(f"  a2={mp.nstr(sol[0],30)}")
print(f"  a4={mp.nstr(sol[1],30)}")
print(f"  a5={mp.nstr(sol[2],30)}")
# is a2/a4 nice? 8/3=2.6667
print("a2 vs 8/3:", mp.nstr(sol[0]-mp.mpf(8)/3,6))
# residuals
for nm,e in [('C4',C4b),('Q7',Q7b),('Q11',Qb)]:
    v=e.subs({a2:sp.Float(mp.nstr(sol[0],40)),a4:sp.Float(mp.nstr(sol[1],40)),a5:sp.Float(mp.nstr(sol[2],40))})
    print(f"  {nm} resid={mp.nstr(abs(mp.mpf(str(v))),4)}")
