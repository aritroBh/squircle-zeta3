# step_22.py : COMPREHENSIVE verification: C-expansion at e^2 level, R, R3, R4 assembly and residual checks
import math, pickle
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta
mp.dps = 70
L2 = log(2); z2 = pi**2/6; z3 = zeta(3); z4 = pi**4/90; z5 = zeta(5); z6 = pi**6/945

# e_k from E_g = Γ((1+z)/2)Γ(1-z/2)/√π : log E_g = -L z + Σ_{k≥2} c_k z^k, c_k = [(-1)^k(2^k-1)+1]ζ(k)/(k 2^k)
cs = {2: z2/2, 3: -z3/4, 4: z4/4, 5: -3*z5/16, 6: z6/6}
def build_ek(K):
    import sympy as sp
    zz = sp.symbols('zz')
    s = -sp.log(2)*zz + sum(sp.zeta(k)*((-1)**k*(2**k-1)+1)*zz**k/(k*2**k) for k in range(2, K+1))
    ser = sp.series(sp.exp(s), zz, 0, K+1).removeO().expand()
    return [mpf(float(sp.N(ser.coeff(zz,k), 40))) for k in range(K+1)]
eks = build_ek(12)
mt = [sum(eks[:k+1]) for k in range(len(eks))]   # m̃_k
def poly(k, N):   # poly_k(g)(N) = L^{k+1}/(k+1) + Σ_{i=1}^{k} k!/(k+1-i)! m̃_i L^{k+1-i} + k! m̃_{k+1}
    L = log(mpf(N)); s = L**(k+1)/mpf(k+1)
    for i in range(1, k+2):
        s += mpf(math.factorial(k))/mpf(math.factorial(k+1-i))*mt[i]*L**(k+1-i)
    return s

d = pickle.load(open("E_values_new.pkl","rb"))
print("=== C(n) expansion at O(e^2): C = poly1 + e(poly1-poly2) - e/(2√2) + e²[(3/4)poly1-(5/4)poly2+(1/2)poly3-3/(2√2)] ===")
for k,v in d.items():
    n = mpf(k); C = mpf(v[2]); e = 1/n; N = n*sqrt(2)
    p1,p2,p3 = poly(1,N),poly(2,N),poly(3,N)
    Cpred = p1 + e*(p1-p2) - e/(2*sqrt(2)) + e**2*(mpf('0.75')*p1 - mpf('1.25')*p2 + mpf('0.5')*p3 - 3/(2*sqrt(2)))
    print("n=%s: C=%s  Cpred=%s  diff=%s" % (k, nstr(C,32), nstr(Cpred,32), nstr(C-Cpred,8)))

# R4 assembly: R4(n) = -3√2/4 - (3/4)poly1 + (5/4)poly2 - (1/2)poly3  (as poly in L_N)
def R4pred(n):
    N = mpf(n)*sqrt(2); p1,p2,p3 = poly(1,N),poly(2,N),poly(3,N)
    return -3*sqrt(2)/4 - mpf('0.75')*p1 + mpf('1.25')*p2 - mpf('0.5')*p3
# R3 with prev d3
b3 = mpf('0.5') - L2/2
c3 = L2**2/4 - L2/2 + 1 + z2
d3 = -L2**3/24 + L2**2/8 - L2/2 - L2*z2/2 + z2/2 - z3/2 + 1 + 3*sqrt(2)/2
def R3(n):
    L = log(mpf(n))
    return L**3/3 + b3*L**2 + c3*L + d3
beta = L2/2 - 1
gamma_ = -(1+sqrt(2)) + L2/2 - pi**2/12 - L2**2/8
print()
print("=== E-residual chain: resid = n²(E-√2/n)-R(n) ~ R3(n)/n ; n(resid-R3(n)/n) = R3num-R3pred ~ R4(n)/n ===")
for k,v in d.items():
    n = mpf(k); E = mpf(v[0]); resid = mpf(v[1])
    L = log(n)
    ratio = resid/(L**3/n)
    R3num = n*resid; d3diff = R3num - R3(n)
    R4 = d3diff*n
    print("n=%s: resid=%s ratio[(log^3)/n]=%s" % (k, nstr(resid,22), nstr(ratio,10)))
    print("   R3num-R3pred=%s  R4pred=%s  diff=%s" % (nstr(d3diff,14), nstr(R4pred(n),14), nstr(R4-R4pred(n),8)))
