import sys, importlib
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp
from mk_integrands import t_integrand
import genintegrand; importlib.reload(genintegrand)
from genintegrand import S_poly, w, binom_half

t = sp.Symbol('t', positive=True)
q2 = t**2 + (1-t)**2
q = sp.sqrt(q2)
ti = t_integrand(2)
pref = -sp.Rational(1,2)/sp.factorial(2)
cM = sp.simplify(pref*ti[(1,1)])
print("ref (1,1):", cM)
S = S_poly(2,1)
print("S_1(t) =", S)
Ssym = sum(sp.Rational(s_) * t**n for n, s_ in enumerate(S))
mine = q*Ssym/(2*t*(1-t)*q2**2)
print("mine:", sp.simplify(mine))
print("ratio:", sp.simplify(sp.together(cM/mine)))
# Check the u-expansion numerically: coefficient of eps^2 l1 l2 in -sqrt(1+u)
# u = alpha X + beta Y; sqrt(1+u) = 1 + u/2 - u^2/8 + ...
# [eps^2 l1 l2] u/2 = [eps^2 l1 l2] (alpha X + beta Y)/2 = 0  (X,Y have pure powers)
# [eps^2 l1 l2] (-u^2/8) = -2 alpha beta [eps l1]X [eps l2]Y /8 = -2 alpha beta (2)(2)/8 = -alpha beta
# so [eps^2]Phi has l1 l2 coeff = -Q'*(-alpha beta) = Q' alpha beta;  M-coeff = -Q' alpha beta/2
ab = (1-t)**2/q2 * t**2/q2
print("expected:", sp.simplify(-q/(t*(1-t))/2 * ab))
print("Q' alpha beta /2 vs S-formula:", sp.simplify(mine - (-q/(2*t*(1-t))*ab)))
