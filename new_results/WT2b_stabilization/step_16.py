import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, ".")
import sympy as sp

v, eps, l1, l2 = sp.symbols('v eps l1 l2')
S2 = sp.sqrt(2)
t = (v**2+2*v-1)/(4*v)
q2 = t**2 + (1-t)**2
q = (v**2+1)/(2*S2*v)
print("q2 == q^2 as rational fns:", sp.simplify(q2 - q**2) == 0)

# alpha + beta == 1 ?
alpha = (1-t)**2/q2; beta = t**2/q2
print("alpha+beta == 1:", sp.simplify(alpha+beta) == 0)

# The q-block of the M-integrand generating series:
#   Qq(eps;v) = (q/(2 t (1-t))) sqrt(1+u) dt/dv,  1+u = alpha e^{2 eps l1} + beta e^{2 eps l2}
# Claim: equals sqrt((1-t)^2 e^{2 eps l1} + t^2 e^{2 eps l2})/(2 t (1-t)) * dt/dv
u = alpha*(sp.exp(2*eps*l1)-1) + beta*(sp.exp(2*eps*l2)-1)
one_plus_u = sp.simplify(1+u)
print("1+u == ( (1-t)^2 e^{2el1} + t^2 e^{2el2} )/q2 :",
      sp.simplify(one_plus_u - ((1-t)**2*sp.exp(2*eps*l1)+t**2*sp.exp(2*eps*l2))/q2)==0)
# hence q*sqrt(1+u) = sqrt(q2)*sqrt(1+u) = sqrt(q2(1+u)) = sqrt((1-t)^2 e^{2el1}+t^2 e^{2el2})
# (branch: as eps-series both sides have constant term 1 and q=sqrt(q2) is the rational determination)
print("=> q-block 1-form:  omega_q = sqrt((1-t)^2 e^{2eps l1}+t^2 e^{2eps l2})/(2t(1-t)) dt")
# verify numerically at a generic point against series expansion
import mpmath as mp
mp.mp.dps = 30
tv = (mp.mpf('0.7')**2+2*mp.mpf('0.7')-1)/(4*mp.mpf('0.7'))
epv = mp.mpf('0.05'); L1 = mp.log(tv); L2 = mp.log(1-tv)
q2v = tv**2+(1-tv)**2; qv = mp.sqrt(q2v)
al = (1-tv)**2/q2v; be = tv**2/q2v
lhs = qv*mp.sqrt(1+al*(mp.e**(2*epv*L1)-1)+be*(mp.e**(2*epv*L2)-1))
rhs = mp.sqrt((1-tv)**2*mp.e**(2*epv*L1)+tv**2*mp.e**(2*epv*L2))
print("numeric check q*sqrt(1+u) vs sqrt(...):", mp.nstr(lhs-rhs,5))
