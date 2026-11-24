import sys, os
workdir = '<project>/research/08_theory/TOW2R_review'
sys.path.insert(0, os.path.join(workdir, 'tow2_copy'))
os.chdir(os.path.join(workdir, 'tow2_copy'))
import mpmath as mp
mp.mp.dps = 50
from mpmath import sqrt, exp, pi

def t_of_v(v): return (v**2+2*v-1)/(4*v)
def q2_of_t(t): return t**2+(1-t)**2

def Phi_eps(eps, t, l1, l2):
    q2 = q2_of_t(t); q = sqrt(q2)
    X = exp(2*eps*l1)-1; Y = exp(2*eps*l2)-1
    alpha = (1-t)**2/q2; beta = t**2/q2
    u = alpha*X+beta*Y
    Qp = q/(t*(1-t))
    A = 1; B = 1   # A,B blocks are constants in eps expansion of exp blocks; we only need sqrt part
    return -Qp*sqrt(1+u)

# sqrt part of -1/2 [eps^d] Phi, coefficient of l1^m l2^(d-m), as f(t); then * dt/dv at v
def sqrt_coeff_numeric(v, d, m):
    t = t_of_v(v)
    ser = mp.taylor(lambda e: Phi_eps(e, t, 1.0, 1.0), 0, d)  # in eps with l1=l2=1: NOT separating l's
    return ser

# Better: coefficient of l1^m l2^(d-m) in [eps^d](-1/2)(-Qp sqrt(1+u)):
# u = eps*(2*(alpha l1 + beta l2)) + eps^2*(2*(alpha l1^2 + beta l2^2)) + ...
# sqrt(1+u) = sum_j binom(1/2,j) u^j ; u^j = sum over compositions.
# Direct numeric: use sympy-free series in 3 vars? Use mpmath with dummy vars via taylor in eps only, with l1, l2 formal:
import itertools
def u_series_coeffs(t, d):
    # returns U[n][(a,b)] = coeff of eps^n l1^a l2^b in u
    q2 = q2_of_t(t); alpha=(1-t)**2/q2; beta=t**2/q2
    U = {}
    for n in range(1, d+1):
        # X^s Y^{j-s}: X = sum_{r>=1} (2 eps l1)^r/r!
        pass
    return U

# Simplest: multinomial expansion. u^j = (alpha X + beta Y)^j = sum_s binom(j,s) alpha^s beta^{j-s} X^s Y^{j-s}
# X^s = sum_{r>=s} S(r,s) s! (2 eps l1)^r/r! ... using (e^z-1)^s = s! sum_{r>=s} S(r,s) z^r/r!
from math import comb, factorial
def S2(n,k):
    if n<0 or k<0: return mp.mpf(0)
    if n==0 and k==0: return mp.mpf(1)
    if k==0 or n==0 or k>n: return mp.mpf(0)
    return k*S2(n-1,k)+S2(n-1,k-1)

def M_k_sqrt_at_v(v, d, m):
    """the sqrt part coefficient of l1^m l2^(d-m) in (-1/2)[eps^d]Phi, times dt/dv, at v."""
    t = t_of_v(v); q2 = q2_of_t(t); q = mp.sqrt(q2)
    alpha=(1-t)**2/q2; beta=t**2/q2
    tot = mp.mpf(0)
    for j in range(1, d+1):
        bj = mp.binomial(mp.mpf('0.5'), j)
        for s in range(0, j+1):
            # coeff of l1^m l2^(d-m) eps^d in (alpha X)^s (beta Y)^{j-s}:
            # X^s contributes s! S(m,s) (2eps)^m/m! ... need m >= s; Y^{j-s}: (j-s)! S(d-m, j-s)(2eps)^{d-m}/(d-m)!
            if m < s or d-m < j-s: continue
            cx = factorial(s)*S2(m,s)*2**m/factorial(m)
            cy = factorial(j-s)*S2(d-m,j-s)*2**(d-m)/factorial(d-m)
            tot += bj*comb(j,s)*alpha**s*beta**(j-s)*cx*cy
    Qp = q/(t*(1-t))
    dtdv = (v**2+1)/(4*v**2)
    return mp.mpf('0.5')*Qp*tot*dtdv   # (-1/2)*(-Qp)*... 

# compare with -q S(t)/(2 t(1-t) q2^d) * dt/dv using genintegrand S_poly
import genintegrand as GI
from fractions import Fraction as Fr
def Sval(S, t):
    return sum(mp.mpf(c.numerator)/c.denominator * t**k for k,c in enumerate(S))
def formula_coeff(v, d, m):
    t = t_of_v(v); q2 = q2_of_t(t); q = mp.sqrt(q2)
    S = GI.S_poly(d, m)
    dtdv = (v**2+1)/(4*v**2)
    return -q*Sval(S,t)/(2*t*(1-t)*q2**d)*dtdv

v = mp.mpf('0.7')+mp.mpf('0.13')*1j
for (d,m) in [(4,2),(4,4),(5,3)]:
    a = M_k_sqrt_at_v(v,d,m); b = formula_coeff(v,d,m)
    print((d,m), 'match:', abs(a-b) < mp.mpf('1e-40'), abs(a-b))
