import sympy as sp
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
lines=open('R_factors.txt').read().split('\n')
G=None
for l in lines:
    if 'deg 22' in l:
        G=sp.sympify(l.split(':',1)[1].strip())
print("G terms:", len(sp.Add.make_args(G)), "variables:", G.free_symbols)
# verify G is in (a1,a2,a4) only (a3=1 dehomogenized)
A_expr=sp.sympify(open('subres1_A.txt').read())
B_expr=sp.sympify(open('subres1_B.txt').read())
print("A vars:", A_expr.free_symbols, " B vars:", B_expr.free_symbols)
# sanity: G(1,1,a4) should have roots ~1.2937 and ~5.1967
import mpmath as mp
mp.mp.dps=50
g1=sp.lambdify(a4, G.subs({a1:1,a2:1}), 'mpmath')
# find positive roots
import numpy as np
Gpoly=sp.Poly(G.subs({a1:1,a2:1}), a4)
print("deg of G(1,1,a4) in a4:", Gpoly.degree())
coeffs=[mp.mpf(int(c)) for c in Gpoly.all_coeffs()]
roots=mp.polyroots(coeffs)
proots=[r for r in roots if abs(sp.im(str(r)))<1e-30 and sp.re(str(r))>0]
print("positive real roots of G(1,1,a4):")
for r in roots:
    if abs(mp.im(r))<mp.mpf('1e-25') and mp.re(r)>0:
        print("   ", mp.nstr(mp.re(r),20))
