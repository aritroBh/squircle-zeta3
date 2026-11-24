import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import t_integrand
import genintegrand; importlib.reload(genintegrand)
from genintegrand import S_poly, Q2

t = sp.Symbol('t', positive=True)
ti = t_integrand(2)
pref = -sp.Rational(1,2)/sp.factorial(2)
q2 = t**2 + (1-t)**2
q = sp.sqrt(q2)
for (i,j), c in sorted(ti.items(), reverse=True):
    cM = sp.simplify(pref * c)   # M_2 coefficient of l1^i l2^j
    S = S_poly(2, i)
    Ssym = sum(sp.Rational(s) * t**n for n, s in enumerate(S))
    mine = -q*Ssym/(2*t*(1-t)*q2**2)
    if i == 2: mine = mine - 1/(2*2*t)
    if i == 0: mine = mine - 1/(2*2*(1-t))
    print((i,j), "t-level match:", sp.simplify(sp.together(cM - mine)) == 0)
    if sp.simplify(sp.together(cM - mine)) != 0:
        print("   ref:", sp.simplify(cM))
        print("   mine:", sp.simplify(mine))
