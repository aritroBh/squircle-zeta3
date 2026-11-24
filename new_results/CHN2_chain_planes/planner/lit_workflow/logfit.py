import sys
sys.path.insert(0,"<project>/research/08_theory/CHN1_chain_bodies")
from fractions import Fraction as F
import sympy as sp, engine, mzv_exact as M
z2,z3,z5,z7,z9,z11=sp.symbols('z2 z3 z5 z7 z9 z11')
even={2:z2,4:sp.Rational(2,5)*z2**2,6:sp.Rational(8,35)*z2**3,8:sp.Rational(24,175)*z2**4,10:sp.Rational(32,385)*z2**5}
zs={2:z2,3:z3,5:z5,7:z7,9:z9,11:z11}
def zval(n): return even[n] if n%2==0 else zs[n]
def basis_to_expr(W,coeffs):
    expr=0; bad=[]
    for mon,c in zip(M.BASIS[W],coeffs):
        if c==0: continue
        if len(mon)==1 and len(mon[0])>1:   # depth>1 irreducible
            bad.append((mon,c)); continue
        term=sp.Integer(1)
        for fac in mon:
            assert len(fac)==1
            term*=zval(fac[0])
        expr+=sp.Rational(c.numerator,c.denominator)*term
    return expr,bad
def logseries(a,Wmax):
    r=engine.reduced_coeffs(tuple(F(x) for x in a),Wmax)
    u=sp.symbols('u')
    ser=1; allbad={}
    for W in range(2,Wmax+1):
        e,bad=basis_to_expr(W,r[W]); ser+=e*u**W
        if bad: allbad[W]=bad
    # the engine series is Vol/16 expansion; coefficient c1 (weight1) is 0 here? include
    L=sp.series(sp.log(ser),u,0,Wmax+1).removeO()
    return sp.expand(L),allbad,u
for a in [(1,1,2,2),(2,2,1,1),(1,2,1,2),(1,3,1,3)]:
    L,bad,u=logseries(a,9)
    print("a=",a,"irreducible leftovers:",bad)
    print("  log coefficients by weight:")
    for W in range(2,10):
        print("   W=%d:"%W, sp.collect(sp.expand(L.coeff(u,W)),[z2,z3,z5,z7,z9]))
