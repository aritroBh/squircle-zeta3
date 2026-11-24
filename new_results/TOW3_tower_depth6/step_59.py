import sys, os, pickle
from fractions import Fraction as Fr
work = '<project>/research/08_theory/TOW3_tower_depth6'
sys.path.insert(0, work); os.chdir(work)
import sympy as sp
import tow3lib as T3
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def g_of_expr(expr, dd, mm):
    v = sp.expand(expr.subs({T3.d: dd, T3.m: mm}))
    return (sp.Rational(sp.re(v)), sp.Rational(sp.im(v)))
dd = 2; delta = 1; s = 0
tot = (Fr(0),Fr(0))
for dpp in range(0, delta+1):
    g = g_of_expr(T3.e(dpp,s), 2*dd+1, dd+1)
    print('dpp',dpp,'g=e-ratio=',g)
    t1 = T3.tau1(delta-dpp)
    t1g = (sp.Rational(sp.re(t1)), sp.Rational(sp.im(t1)))
    print('  tau1 =', t1g)
    rB0d = (Fr(0), Fr(2*(4*dd-1), dd+1))
    term = gmul(gmul(g, t1g), rB0d)
    term = (term[0]*Fr(dd+1,4*dd-dpp-1), term[1]*Fr(dd+1,4*dd-dpp-1))
    print('  term =', term)
    tot = (tot[0]+term[0], tot[1]+term[1])
print('ehat(1,0) =', tot, ' expected (-2,-2)')
