step_03 = r'''# step_03.py -- D11c(1) rhombus: kappa_n(a,b) closed form, symbolic verification
import sympy as sp
from sympy import symbols, Rational, factorial, binomial, sqrt, exp, series

NMAX = 12
Z = {n: symbols(f'Z{n}') for n in range(2, NMAX+1)}
eta, s, t, q, z = symbols('eta s t q z')

def ts_mul(A, B, N=NMAX):
    R = {}
    for i, aa in A.items():
        for j, bb in B.items():
            if i + j <= N: R[i+j] = R.get(i+j, 0) + aa*bb
    return {k: sp.expand(v) for k, v in R.items() if v != 0}
def ts_exp(A, N=NMAX):
    R = {0: sp.Integer(1)}
    for n in range(1, N+1):
        acc = 0
        for k in range(1, n+1):
            if k in A and (n-k) in R: acc += k*A[k]*R[n-k]
        v = sp.expand(acc/n)
        if v != 0: R[n] = v
    return R
def bs_exp(A, N=NMAX):
    R = {(0,0): sp.Integer(1)}
    for tot in range(1, N+1):
        for i in range(0, tot+1):
            j = tot-i
            if i > 0:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 >= 1 and (i-i1, j-j1) in R: acc += i1*aa*R[(i-i1,j-j1)]
                v = sp.expand(acc/i)
            else:
                acc = 0
                for (i1,j1), aa in A.items():
                    if i1 == 0 and j1 >= 1 and (0, j-j1) in R: acc += j1*aa*R[(0,j-j1)]
                v = sp.expand(acc/j) if j > 0 else sp.Integer(1)
            if v != 0: R[(i,j)] = v
    return R

# even part of sqrt(cosh z + q sinh z)
gser = series(sqrt(sp.cosh(z) + q*sp.sinh(z)), z, 0, NMAX+1).removeO()
cq = {k: sp.expand((gser.coeff(z,2*k))) for k in range(NMAX//2+1)}  # odd coeffs irrelevant (killed by symmetry)
cq = {k: sp.expand(cq[k]) for k in cq}
print("c_k(q) =", {k: cq[k] for k in range(4)})

mA = {n: (-1)**(n+1)*Z[n]/n for n in range(2, NMAX+1)}
K = ts_exp(mA)
Bser = {}
for n in range(2, NMAX+1):
    for r in range(0, n+1, 2):
        Bser[(n-r,r)] = Bser.get((n-r,r), 0) + (-1)**n*Z[n]*2*binomial(n,r)/n
Bser = {k: sp.expand(v) for k, v in Bser.items() if v != 0}
Hser = bs_exp(Bser)
Sq = {}
for k in range(NMAX//2+1):
    for i in range(0, NMAX+1):
        v = Hser.get((i,2*k), 0)
        if v != 0: Sq[2*k+i] = Sq.get(2*k+i,0) + cq[k]*factorial(2*k)*v/2**i
Sq = ts_mul(K, Sq)

A_, B_ = (1+q)/2, (1-q)/2
Rq = series((sqrt(A_*exp(4*t)+B_) + sqrt(A_+B_*exp(4*t)))/2, t, 0, NMAX+1).removeO()
print(f"\n{'n':>3} {'kappa_n(q) from S_q':>42} {'closed form':>42} ok")
for n in range(3, 12, 2):
    kq = sp.expand(Sq.get(n,0)).coeff(Z[n],1)*4*(-1)**n
    Vq = sp.expand(factorial(n)*Rq.coeff(t,n))
    kf = sp.expand(Rational(4*(-1)**n,n)*(1-Rational(2)**(1-n)*Vq))
    print(f"{n:>3} {str(kq):>42} {str(kf):>42} {sp.simplify(kq-kf)==0}")

V3 = sp.expand(factorial(3)*Rq.coeff(t,3))
print("\nV_3(q) =", sp.factor(V3), "; kappa_3(q) =", sp.factor(sp.expand(Rational(-4,3)*(1-Rational(1,4)*V3))),
      "= 4 a^2 b^2/(a^2+b^2)^2")
# a_3(q) = (1-q^2)(Z3 - 2 Z2): check from Sq
s2 = sp.expand(Sq.get(2,0)); s3 = sp.expand(Sq.get(3,0))
a3q = sp.expand(4*(-2*s2 - s3))
print("a_3(q) =", a3q, " ; (1-q^2)(Z3-2Z2) =", sp.expand((1-q**2)*(Z[3]-2*Z[2])))
'''
open('step_03.py','w').write(step_03)

step_04 = r'''# step_04.py -- D11c(1) rhombus: numerical check of a_3(q)=(1-q^2)(zeta(3)-pi^2/3) via direct arclength
import mpmath as mp
mp.mp.dps = 40

def perim_rhombus(p, a, b):
    def f(th):
        c, s = mp.cos(th), mp.sin(th)
        dxd = a*(2/p)*c**(2/p-1)*(-s); dyd = b*(2/p)*s**(2/p-1)*c
        return mp.sqrt(dxd**2+dyd**2)
    return 4*mp.quad(f, [0, mp.pi/2])

def a3_numeric(a, b, r):
    est = []
    for h in [mp.mpf('0.02'), mp.mpf('0.01'), mp.mpf('0.005')]:
        f2 = perim_rhombus(1+2*h,a,b); f1 = perim_rhombus(1+h,a,b)
        fm1 = perim_rhombus(1-h,a,b); fm2 = perim_rhombus(1-2*h,a,b)
        est.append((f2-2*f1+2*fm1-fm2)/(2*h**3)/6/r)
    e1 = est[0]+(est[0]-est[1])/3; e2 = est[1]+(est[1]-est[2])/3
    return e2+(e2-e1)/15

for (a,b) in [(1,1),(2,1),(3,2)]:
    r = mp.sqrt(a**2+b**2); q = (a**2-b**2)/(a**2+b**2)
    a3 = a3_numeric(a,b,r)
    pred = (1-q**2)*(mp.zeta(3)-mp.pi**2/3)
    print(f"a={a},b={b}: a3 numeric={mp.nstr(a3,18)}  predicted (1-q^2)(zeta3-pi^2/3)={mp.nstr(pred,18)}")
'''
open('step_04.py','w').write(step_04)
print("step_03.py, step_04.py written")
