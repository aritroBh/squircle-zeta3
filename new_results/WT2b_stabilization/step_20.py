import sys, importlib, pickle, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import sympy as sp
import fastred; importlib.reload(fastred)
from fastred import *
from fastred import v as vv

v, eps, l1, l2 = sp.symbols('v eps l1 l2')
S2 = sp.sqrt(2)
t_ = (v**2+2*v-1)/(4*v)
q_ = (v**2+1)/(2*S2*v)
q2_ = t_**2 + (1-t_)**2
dtv_ = sp.diff(t_, v)

# source m(eps;v) = (1/2) sqrtS * dt/dv,  S = t^{-2} e^{2 eps l1} + (1-t)^{-2} e^{2 eps l2}
# [eps^k] sqrtS = q/(t(1-t)) * sigma_k,  sigma_k = [eps^k] (1+sum c_n eps^n)^{1/2}
# c_n = (2^n/n!) ((1-t)^2 l1^n + t^2 l2^n)/q2
KMAX = 8
c = {}
for n in range(1, KMAX+1):
    c[n] = sp.Rational(2**n, sp.factorial(n)) * ((1-t_)**2 * l1**n + t_**2 * l2**n)/q2_
# sigma_k via binomial series
ser = sum(c[n]*eps**n for n in range(1, KMAX+1))
sq = sp.series(sp.sqrt(1+ser), eps, 0, KMAX+1).removeO().expand()
sigma = {k: sp.expand(sq).coeff(eps, k) for k in range(KMAX+1)}
# m_k = (1/2)(q/(t(1-t))) dt/dv sigma_k ; coefficient of l1^a l2^b
pref = sp.Rational(1,2) * q_/(t_*(1-t_)) * dtv_
m = {}
for k in range(KMAX+1):
    poly = sp.Poly(sp.expand(sigma[k]), l1, l2)
    for (a,b), cf in poly.terms():
        m[(k,a,b)] = sp.together(pref*cf)
print("source coefficients m_k^{ab} built for k<=", KMAX, "; count:", len(m))
# sanity: m_0 = (1/2) q/(t(1-t)) dt/dv
print("m_0 has +-i poles? denom factors:", sp.factor(sp.denom(m[(0,0,0)])))
# m_1 should have only SIMPLE +-i poles
d1 = sp.denom(sp.together(m[(1,1,0)]))
print("m_1^10 denom:", sp.factor(d1))
pickle.dump({'KMAX':KMAX}, open('tmp_src_meta.pkl','wb'))
# store m as sympy strings? keep in memory for next step within same call - but calls are fresh processes!
# -> must pickle. pickle sympy exprs via srepr
pickle.dump({k: sp.srepr(val) for k,val in m.items()}, open('tmp_src.pkl','wb'))
print("pickled source")
