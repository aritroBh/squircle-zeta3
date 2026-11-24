"""v_proofs_sympy.py -- symbolic checks for two short proofs written in verify.md (verifier's own).

(P-W2) W2(c) = -int_1^{t1} log g dt/(t+1), g = (t^2-1)/(2(t-ca)), ca = 2c^2-1, t1 = 1+2s, s = sqrt(1-c^2).
   log g(t1) = 0, so dW2/dc = -int_1^{t1} d_c log g dt/(t+1) = -int 4c/((t-ca)(t+1)) dt.
   Check: this equals d/dc[ z2/2 - Li2(c^2)/2 ] = log(1-c^2)/c, and W2(c=0) = pi^2/12 (numerically).
(P-B1) all-n proof that the HPL form of B1 has Taylor coefficients mu_{1,n} = D_n/(n(2n-1)) - 1/(2n^2):
   with A_m = sum_{k<=m} (-1)^(k+1)/k (A_{2n-1} = D_n), the c^N coefficient is -X_N + X_{N-1} + Y_N + Y_{N-1},
   X_N = 1/N^2 + A_{N-1}/N,  Y_N = (-1)^(N+1)/N^2 - B_{N-1}/N,  B_{N-1} = sum_{k<N} (-1)^(N-1-k)/k
   (B_{2n-1} = A_{2n-1}, B_{2n-2} = -A_{2n-2}, B_{2n} = -A_{2n}).
"""
import sympy as sp
import mpmath as mp

c, t, n, A, A2 = sp.symbols('c t n A A2', positive=True)
ca = 2*c**2 - 1
s = sp.sqrt(1 - c**2)
t1 = 1 + 2*s
g = (t**2 - 1)/(2*(t - ca))
print('[P-W2] log g(t1) = 0 :', sp.simplify(g.subs(t, t1) - 1) == 0)
print('       d_c log g - 4c/(t-ca) =', sp.simplify(sp.diff(sp.log(g), c) - 4*c/(t - ca)))
anti = -(4*c/(1 + ca))*(sp.log(t - ca) - sp.log(t + 1))          # antiderivative of -4c/((t-ca)(t+1))
dW2 = sp.simplify(anti.subs(t, t1) - anti.subs(t, 1))
target = sp.log(1 - c**2)/c
print('       dW2/dc - log(1-c^2)/c at c = 0.1, 0.37, 0.8, 0.99:',
      [sp.N(dW2.subs(c, v) - target.subs(c, v), 30) for v in (sp.Rational(1, 10), sp.Rational(37, 100), sp.Rational(4, 5), sp.Rational(99, 100))])
print('       symbolic simplify(dW2 - target) =', sp.simplify(sp.expand_log(dW2 - target, force=True)))
mp.mp.dps = 40
W0 = -mp.quad(lambda tt: mp.log((tt - 1)/2)/(tt + 1), [1, 3])
print('       W2(c=0) - pi^2/12 =', mp.nstr(W0 - mp.pi**2/12, 3))

# P-B1: even N = 2n, A = A_{2n-1}, A_{2n-2} = A - 1/(2n-1)
N = 2*n
X_N = 1/N**2 + A/N
X_Nm1 = 1/(N - 1)**2 + (A - 1/(2*n - 1))/(N - 1)
Y_N = -1/N**2 - A/N                       # B_{2n-1} = A_{2n-1}
Y_Nm1 = 1/(N - 1)**2 + (A - 1/(2*n - 1))/(N - 1)    # -B_{2n-2}/(2n-1) = +A_{2n-2}/(2n-1)
coef_even = sp.simplify(-X_N + X_Nm1 + Y_N + Y_Nm1 - (A/(n*(2*n - 1)) - 1/(2*n**2)))
print('[P-B1] even coefficient minus mu_{1,n} formula:', coef_even)
# odd N = 2n+1, A2 = A_{2n}, A = A_{2n-1}
X_N = 1/(2*n + 1)**2 + A2/(2*n + 1)
X_Nm1 = 1/(2*n)**2 + A/(2*n)
Y_N = 1/(2*n + 1)**2 + A2/(2*n + 1)       # B_{2n} = -A_{2n}
Y_Nm1 = -1/(2*n)**2 - A/(2*n)             # B_{2n-1} = A_{2n-1}
print('       odd coefficient:', sp.simplify(-X_N + X_Nm1 + Y_N + Y_Nm1))
# numeric spot check of the alternating-harmonic facts
from fractions import Fraction as Fr
ok = True
for nv in range(1, 60):
    Aa = lambda m: sum((Fr((-1)**(k + 1), k) for k in range(1, m + 1)), Fr(0))
    Bb = lambda m, Nv: sum((Fr((-1)**(Nv - 1 - k), k) for k in range(1, m + 1)), Fr(0))
    H = lambda m: sum((Fr(1, k) for k in range(1, m + 1)), Fr(0))
    ok &= Aa(2*nv - 1) == H(2*nv - 1) - H(nv - 1)
    ok &= Bb(2*nv - 1, 2*nv) == Aa(2*nv - 1) and Bb(2*nv - 2, 2*nv - 1) == -Aa(2*nv - 2) and Bb(2*nv, 2*nv + 1) == -Aa(2*nv)
print('       A_{2n-1} = H_{2n-1}-H_{n-1} and the B/A sign rules, n < 60:', ok)
