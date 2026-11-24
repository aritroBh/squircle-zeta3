import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import t_integrand
import genintegrand; importlib.reload(genintegrand)
from genintegrand import S_poly

t = sp.Symbol('t', positive=True)
q2 = t**2 + (1-t)**2
q = sp.sqrt(q2)
ok = True
for k in [2,3,4,5]:
    ti = t_integrand(k)
    pref = -sp.Rational(1,2)/sp.factorial(k)
    for (i,j), c in sorted(ti.items(), reverse=True):
        cM = pref * c
        S = S_poly(k, i)
        Ssym = sum(sp.Rational(s_) * t**n for n, s_ in enumerate(S))
        mine = q*Ssym/(2*t*(1-t)*q2**k)       # corrected sign
        if i == k: mine = mine - 1/(2*sp.factorial(k)*t)
        if i == 0: mine = mine - 1/(2*sp.factorial(k)*(1-t))
        if sp.simplify(sp.together(cM - mine)) != 0:
            ok = False; print("MISMATCH", k, (i,j))
print("t-level closed form matches t_integrand for k=2..5:", ok)
