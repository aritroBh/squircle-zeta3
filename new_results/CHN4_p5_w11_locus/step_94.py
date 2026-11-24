import sympy as sp, mpmath as mp, pickle
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
def residuals(a1v,a2v,a4v,a5v):
    dd={a1:sp.Float(mp.nstr(a1v,50)),a2:sp.Float(mp.nstr(a2v,50)),a3:1,a4:sp.Float(mp.nstr(a4v,50)),a5:sp.Float(mp.nstr(a5v,50))}
    return max(abs(mp.mpf(str(e.subs(dd)))) for e in (C4,Q7,Q11sp))
def gamma_point(a1v, guess):
    C4s=C4.subs({a1:sp.Rational(a1v),a3:1}); Q7s=Q7.subs({a1:sp.Rational(a1v),a3:1}); Qs=Q11sp.subs({a1:sp.Rational(a1v),a3:1})
    f=[sp.lambdify((a2,a4,a5), e, 'mpmath') for e in (C4s,Q7s,Qs)]
    return mp.findroot(f, guess, tol=mp.mpf('1e-45'), maxsteps=400)
# continue down from a1=0.0005 (down branch ended at a2=2.67043,a4=2.67045,a5=1.00233)
down=pickle.load(open('gamma_down.pkl','rb'))
cur=(down[-1][1],down[-1][2],down[-1][3])
newsmall=[]
for a1v in [0.0002, 0.0001, 0.00005]:
    try:
        sol=gamma_point(sp.Rational(a1v), cur); 
        if all(x>0 for x in sol):
            cur=(sol[0],sol[1],sol[2]); newsmall.append((mp.mpf(a1v),cur[0],cur[1],cur[2]))
    except Exception as ex:
        print(f"fail a1={a1v}"); break
for r in newsmall:
    print(f"  a1={mp.nstr(r[0],7)} a2={mp.nstr(r[1],8)} a4={mp.nstr(r[2],8)} a5={mp.nstr(r[3],8)} resid={mp.nstr(residuals(*r),3)}")
pickle.dump(newsmall, open('gamma_small.pkl','wb'))
print("saved gamma_small.pkl")
