import sympy as sp, mpmath as mp, pickle
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
x,y=sp.Rational(1,100), sp.Rational(3)
Gp=sp.Poly(G.subs({a1:x,a2:y}), a4, domain=sp.QQ)
ivs=Gp.intervals()
print("isolating intervals for a4 (lo,hi,mult):")
for (lo,hi),m in ivs:
    print(f"  ({float(lo):.5f},{float(hi):.5f}) mult={m} {'POS' if hi>0 else ''}")
