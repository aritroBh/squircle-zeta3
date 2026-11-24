"""V4. Asymptotic expansion of L_theta(n) through n^-4, derived symbolically from the exact formula and
checked numerically against the exact formula (itself checked against quadrature in V2) at large n.

  L_theta(n) = 1 + c + 2z (kappa/n)^(1/(n-1)) + Q2/n^2 + Q3(log n)/n^3 + Q4(log n)/n^4 + O(log^3 n / n^5)
  kappa = e^{a1},  log R_theta(eps) = log(2z) + a1 eps + a2 eps^2 + ...,
  S(n) = sum_m gamma_m eps n^(1-m)/(m-1-eps).
"""
import sympy as sp
import mpmath as mp
from common import *
from closed_forms import closed_a

# ---------- symbolic expansion ----------
N, Lg, z = sp.symbols('N Lg z')   # N = 1/n, Lg = log n
a1, a2, a3, a4, a5 = sp.symbols('a1:6')
g2, g3, g4, g5 = sp.symbols('g2:6')
K = 6
eps = N/(1 - N)
E = sp.exp(a2*eps**2 + a3*eps**3 + a4*eps**4 + a5*eps**5)
expo = sp.exp(-eps*(Lg - a1))
T = 2*z*expo*(E - 1)
S = sum(g*eps*N**(m - 1)/(m - 1 - eps) for m, g in zip(range(2, 6), [g2, g3, g4, g5]))
ser = sp.expand(sp.series(T - S, N, 0, K).removeO())
Q = {k: sp.collect(sp.expand(ser.coeff(N, k)), Lg) for k in range(2, K)}
print("Symbolic coefficients of n^-k in L_theta(n) - [1 + c + 2z (e^a1/n)^(1/(n-1))]   (Lg = log n, g_m = gamma_m):")
for k in range(2, 5):
    print(f"  Q{k} =", Q[k])

# with gamma_m substituted: g2 = 2z(1-z), g3 = -2cz(1-z), g4 = (P2 - P4)/7, c = 1-2z
cc = 1 - 2*z
P = [sp.Integer(1), cc]
for k in range(1, 5):
    P.append(sp.expand(((2*k + 1)*cc*P[k] - k*P[k - 1])/(k + 1)))
gam = {m: sp.factor(sp.expand((-1)**m*(P[m - 2] - P[m])/(2*m - 1))) for m in range(2, 6)}
print("\n  gamma_m in terms of z = sin^2(theta/2):", {m: gam[m] for m in gam})
Qz = {k: sp.collect(sp.expand(Q[k].subs({g2: gam[2], g3: gam[3], g4: gam[4], g5: gam[5]})), [Lg]) for k in Q}
print("\n  Q2 =", sp.factor(Qz[2]))
print("  Q3 =", sp.collect(sp.factor_terms(Qz[3]), Lg))
print("  Q4 =", sp.collect(sp.expand(Qz[4]), Lg))

# ---------- numerical check ----------
setdps(100)
fQ = {k: sp.lambdify((Lg, z, a1, a2, a3, a4, a5), Qz[k], 'mpmath') for k in Qz}
angles = [('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3), ('3pi/4', 3*mp.pi/4),
          ('pi/6', mp.pi/6), ('5pi/6', 5*mp.pi/6)]
print("\nNumerical check (dps=100). resid_k = L_exact - [main + sum_{j<=k} Q_j/n^j];"
      " expected sizes ~ log^(k-1)n / n^(k+1)")
for name, th in angles:
    zz = mp.sin(th/2)**2
    cth = mp.cos(th)
    A = closed_a(th)
    a5num = mp.mpf(0)   # a5 enters only Q5 and beyond
    print(f" theta={name}:  a1..a4 = {[mp.nstr(x, 12) for x in A]}")
    for n in [mp.mpf(10)**6, mp.mpf(10)**10, mp.mpf(10)**14]:
        Lx = L_exact(n, th)
        lg = mp.log(n)
        main = 1 + cth + 2*zz*mp.power(mp.exp(A[0])/n, 1/(n - 1))
        r = Lx - main
        out = []
        for k in range(2, 5):
            r -= fQ[k](lg, zz, *A, a5num)/n**k
            out.append(r)
        print(f"   n=1e{int(mp.log10(n)):<2d}  L={mp.nstr(Lx,35)}  resid2={mp.nstr(out[0],4):>10}  resid3={mp.nstr(out[1],4):>10}"
              f"  resid4={mp.nstr(out[2],4):>10}   resid4*n^5/log^3 n={mp.nstr(out[2]*n**5/lg**3,5)}")
    # empirical Q-values at n = 1e14 (reported for the table)
    n = mp.mpf(10)**14
    lg = mp.log(n)
    print(f"   Q2 = {mp.nstr(fQ[2](lg, zz, *A, 0), 40)}")
