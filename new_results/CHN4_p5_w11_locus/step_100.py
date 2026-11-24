import sympy as sp, mpmath as mp, pickle, os
mp.mp.dps=70
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
C4b=C4.subs({a1:0,a3:1}); Q7b=Q7.subs({a1:0,a3:1}); Qb=Q11sp.subs({a1:0,a3:1})
f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4b,Q7b,Qb)]
sol=mp.findroot(f, (mp.mpf('2.667'),mp.mpf('2.667'),mp.mpf('1.0')), tol=mp.mpf('1e-55'), maxsteps=800)
print("boundary point refined: a2=",mp.nstr(sol[0],45))
print("                         a4=",mp.nstr(sol[1],45))
print("                         a5=",mp.nstr(sol[2],45))
print("a2-a4 =",mp.nstr(sol[0]-sol[1],6))
print("\nfinal batch log so far:")
print(open('gamma_final.log').read() if os.path.exists('gamma_final.log') else "not yet")
