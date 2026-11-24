import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import v_integrand
import fastred; importlib.reload(fastred)
from fastred import *

v = sp.Symbol('v', positive=True)
def numeq(e1, e2):
    d = sp.together(e1) - sp.together(e2)
    vals = [complex(sp.N(d.subs(v, x), 25)) for x in (sp.Rational(7,10), sp.Rational(19,10), sp.Rational(13,10))]
    return all(abs(z) < 1e-18 for z in vals)

R = v_integrand(2)[(2,0)]
P = pf_from_sympy(R)
simple, higher, const = pf_split(P)
hs = pf_to_sympy([{}, higher[1], higher[2], ZERO])
# sympy integrate of higher
G_sym = sp.together(sp.integrate(sp.together(hs), v))
G = pf_integrate_higher(higher)
print("integrate ok:", numeq(pf_to_sympy(G), G_sym))
C = pf_eval_at(G, BQ)
print("G(b) mine:", C.sympy(), " sympy:", sp.radsimp(G_sym.subs(v, b)), " equal:", sp.simplify(C.sympy() - G_sym.subs(v,b))==0)
Gm = [G[0], G[1], G[2], -C]
Gm_sym = G_sym - G_sym.subs(v, b)
DL1s = 1/(v-b) + 1/(v+a) - 1/v
prod = pf_mul_DL(Gm, DL1)
print("mul_DL ok:", numeq(pf_to_sympy(prod), sp.together(Gm_sym*DL1s)))
folded = pf_fold(prod)
print("fold ok:", numeq(pf_to_sympy(folded), sp.together((Gm_sym*DL1s).subs(v, 1/v)/v**2)))
