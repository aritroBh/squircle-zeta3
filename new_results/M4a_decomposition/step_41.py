import sys, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
sys.path.insert(0, '<project>/research/08_theory/A3_W4_proof')
import mpmath as mp
mp.mp.dps = 200

# 1) thm:W4 closed form via w4_symbol.conjecture_terms (exact transcription of the theorem)
from w4_symbol import conjecture_terms
from qsqrt2 import QS2, A, B, SQ2
r2 = mp.sqrt(2)
def qnum(c):  # QS2 -> mpf
    return mp.mpf(c.p.numerator)/c.p.denominator + r2*mp.mpf(c.q.numerator)/c.q.denominator
def atomval(a):
    if a[0]=='li':
        z = qnum(a[2])
        v = mp.polylog(a[1], z)
        assert abs(mp.im(v)) < mp.mpf(10)**(-150), (a, v)
        return mp.re(v)
    if a[0]=='log':
        return mp.log(qnum(a[1]))
    if a[0]=='pi2':
        return mp.pi**2
    if a[0]=='z3':
        return mp.zeta(3)
    raise ValueError(a)
W4th = mp.mpf(0)
for c, atoms in conjecture_terms():
    t = qnum(c)
    for a in atoms:
        t *= atomval(a)
    W4th += t
print('W4 theorem value =', mp.nstr(W4th, 120))

# 2) M3 by method-A quadrature (full integrand) at 200 dps
import sympy as sp
from mk_integrands import v_integrand, v as vs
aN = 1+r2; bN = r2-1
Rk = v_integrand(3)
fs = [(i,j,sp.lambdify(vs,R,modules='mpmath')) for (i,j),R in Rk.items()]
def fM3(x):
    l1 = mp.log((x-bN)*(x+aN)/(4*x)); l2 = mp.log((aN-x)*(x+bN)/(4*x))
    return mp.fsum(R(x)*l1**i*l2**j for i,j,R in fs)
M3,err = mp.quad(fM3, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
print('M3 quad =', mp.nstr(M3, 120), 'err', mp.nstr(err,3))

lam = mp.log(2)
T = (mp.pi**2+6*lam**2)/24

# 3) W3' reduced pickle value (k=3)
t3 = pickle.load(open('m3_tails_exact.pkl','rb'))
W3p_red = t3['W3p'].num()
print('W4via reduced  - thm =', mp.nstr(M3 - W3p_red - T - W4th, 6))

# 4) W3' by my quadrature of the k=3 pieces
import json
def quad_pieces(K):
    from ibp_reduce import reduce_k
    import sympy as s2
    out = reduce_k(K, verbose=False)
    a_, b_ = 1+s2.sqrt(2), s2.sqrt(2)-1
    tot = mp.mpf(0)
    for (i,j), R in out.items():
        if i+j != 2: continue
        f = s2.lambdify(vs, R, modules='mpmath')
        g = lambda x: f(x)*mp.log((x-bN)*(x+aN)/(4*x))**i*mp.log((aN-x)*(x+bN)/(4*x))**j
        val,e = mp.quad(g, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
        tot += val
    return tot
W3p_quad = quad_pieces(3)
print('W3p quad =', mp.nstr(W3p_quad, 120))
print('W4via quad     - thm =', mp.nstr(M3 - W3p_quad - T - W4th, 6))
print('reduced - quad (k3)  =', mp.nstr(W3p_red - W3p_quad, 6))
