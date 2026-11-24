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
up2=pickle.load(open('gamma_up2.pkl','rb'))
# refine points near a1=1 and a1=2
targets=[]
for want in [1.0, 2.0]:
    r=min(up2, key=lambda r: abs(float(r[0])-want))
    sol=gamma_point(sp.Rational(float(r[0])), (r[1],r[2],r[3]))
    targets.append((sol and mp.mpf(float(r[0])), sol[0],sol[1],sol[2]))
# a different region: a1=0.123 near the both-signs cell (a1=0.123,a2=1.474); sheets a4=1.123(q<0),a4=0.384(q>0)
# Gamma between them: start guess a2~1.47,a4~0.7,a5~1.0
try:
    s2=gamma_point(sp.Rational(123,1000),(mp.mpf('1.47'),mp.mpf('0.7'),mp.mpf('1.0')))
    targets.append((mp.mpf('0.123'),s2[0],s2[1],s2[2]))
except Exception as ex:
    print("region2 failed", ex)
print("additional Gamma points:")
for r in targets:
    print(f"  a1={mp.nstr(r[0],6)} a2={mp.nstr(r[1],8)} a4={mp.nstr(r[2],8)} a5={mp.nstr(r[3],8)} resid={mp.nstr(residuals(*r),3)}")
pickle.dump(targets, open('gamma_extra.pkl','wb'))
print("saved gamma_extra.pkl")
