code = r'''# step_03.py: Beta-moment table J(n,k,a,b) = int_0^1 t^n (1-t)^k log^a t log^b(1-t) dt, a+b<=5,
# exact in basis {1, z2, z3, z2^2, z2z3, z5}  (zeta(4) = 2 z2^2/5 used).
# Verifies the weight-5 rules: [z5]J = 24*B(n+1,k+1) for a+b=5, a,b>=1, else 0;
# [z2z3]J = -12*B(n+1,k+1) for (a,b) in {(3,2),(2,3)}, else 0.  Numerical check at 40 digits.
import sympy as sp
from functools import lru_cache
u, v = sp.symbols('u v')
z2, z3, z5, g = sp.symbols('z2 z3 z5 g')
z4 = sp.Rational(2, 5)*z2**2   # zeta(4) = pi^4/90 = (2/5) zeta(2)^2

def H(n, r): return sum(sp.Rational(1, j**r) for j in range(1, n))
def lg_series(n):  # lnGamma(n+w)-lnGamma(n), coefficients of w^1..w^5
    return [-g + H(n, 1), (z2 - H(n, 2))/2, (-2*z3 + 2*H(n, 3))/6,
            (6*z4 - 6*H(n, 4))/24, (-24*z5 + 24*H(n, 5))/120]

@lru_cache(maxsize=None)
def J_table(n, k):
    x, y = n + 1, k + 1
    a, b, c = lg_series(x), lg_series(y), lg_series(x + y)
    S = sp.expand(sum(a[j]*u**(j+1) + b[j]*v**(j+1) - c[j]*(u+v)**(j+1) for j in range(5)))
    E = sp.Integer(1); P = sp.Integer(1)
    for j in range(1, 6):
        P = sp.expand(P*S)
        E = sp.expand(E + P/sp.factorial(j))
    B0 = sp.factorial(n)*sp.factorial(k)/sp.factorial(n + k + 1)
    tab = {}
    for A in range(6):
        for Bb in range(6 - A):
            co = sp.expand(E.coeff(u, A).coeff(v, Bb))
            tab[(A, Bb)] = sp.expand(B0*sp.factorial(A)*sp.factorial(Bb)*co)
    return tab
def J(n, k, a, b): return J_table(n, k)[(a, b)]

if __name__ == '__main__':
    import mpmath as mp
    mp.mp.dps = 40
    zv = {z2: sp.pi**2/6, z3: sp.zeta(3), z5: sp.zeta(5)}
    worst = 0; rule5_ok = rule23_ok = True
    for n in range(0, 7):
        for k in range(0, 3):
            for a in range(6):
                for b in range(6 - a):
                    ex = J(n, k, a, b)
                    assert g not in ex.free_symbols, (n, k, a, b, ex)   # Euler gamma cancels
                    val = sp.N(ex.subs(zv), 40)
                    num = mp.quad(lambda t: t**n*(1-t)**k*mp.log(t)**a*mp.log(1-t)**b, [0, 0.5, 1])
                    worst = max(worst, abs(val - num))
                    B0 = sp.Rational(1, (n+k+1)*sp.binomial(n+k, n))
                    c5 = sp.simplify(ex.coeff(z5)); c23 = sp.simplify(ex.coeff(z2*z3))
                    want5 = 24*B0 if (a + b == 5 and a >= 1 and b >= 1) else sp.Integer(0)
                    want23 = -12*B0 if (a, b) in [(3, 2), (2, 3)] else sp.Integer(0)
                    rule5_ok &= (c5 == want5); rule23_ok &= (c23 == want23)
    print("numeric max error over table:", worst)
    print("z5 rule   [z5]J = 24 B(n+1,k+1) iff a+b=5,a,b>=1 else 0 :", rule5_ok)
    print("z2z3 rule [z2z3]J = -12 B(n+1,k+1) iff (a,b)=(3,2)/(2,3) else 0:", rule23_ok)
    print("examples:")
    for key in [(0,0,4,1),(0,0,3,2),(2,0,3,2),(3,1,2,3),(0,0,5,0),(0,0,0,5)]:
        print("  J%s ="%(key,), sp.collect(sp.expand(J(*key)), [z5, z2*z3, z3, z2]))
'''
open('steps/step_03.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_03.py'],capture_output=True,text=True,timeout=3000)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
open('progress.log','a').write("step_01: E5 reduced built (475 terms, deg 5); cross g_ab extracted, poly in t deg 5. step_02: full E5 (1910 terms), g_ab tau/gam-free, matches reduced. step_03 written.\n")
