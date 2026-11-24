import sympy as sp, mpmath as mp, pickle
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
def gamma_point(a1v, guess):
    C4s=C4.subs({a1:sp.Rational(a1v),a3:1}); Q7s=Q7.subs({a1:sp.Rational(a1v),a3:1}); Qs=Q11sp.subs({a1:sp.Rational(a1v),a3:1})
    f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    return mp.findroot(f, guess, tol=mp.mpf('1e-45'), maxsteps=200)
# trace from a1=0.01 (known: a2=2.7413,a4=2.7417,a5=1.0461) upward and downward
g0=(mp.mpf('2.7413'),mp.mpf('2.7417'),mp.mpf('1.0461'))
import numpy as np
# go up
up=[]
cur=g0
for a1v in [0.02,0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.5,2,3,5,7,10,20,50,100]:
    try:
        sol=gamma_point(sp.Rational(a1v), cur); cur=(sol[0],sol[1],sol[2])
        if all(x>0 for x in cur): up.append((mp.mpf(a1v),cur[0],cur[1],cur[2]))
    except Exception as ex:
        print(f"up failed at a1={a1v}: {type(ex).__name__}"); break
print("UP branch (a1,a2,a4,a5):")
for r in up: print(f"  a1={mp.nstr(r[0],5)} a2={mp.nstr(r[1],8)} a4={mp.nstr(r[2],8)} a5={mp.nstr(r[3],8)}")
pickle.dump(up, open('gamma_up.pkl','wb'))
