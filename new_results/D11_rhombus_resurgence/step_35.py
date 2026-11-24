# ======== STEP 3a: reference F values + Pade approximant of Borel series ========
from mpmath import mp, mpf, mpc, quad, sqrt, cosh, exp, log, pi, factorial
import json
mp.dps = 100
Fcoef = [mpf(v) for v in json.load(open('step2_num200.json'))['f200']][:62]

# ---- reference F(eta) by direct quadrature of the l-integral ----
def Lfun(l): return log(2*cosh(l/2))
def F_direct(eta):
    w0 = lambda l: exp(-l)/(1+exp(-l))**2
    f = lambda l: w0(l)*exp(-eta*Lfun(l))*sqrt(cosh(eta*l))
    return quad(f, [-mp.inf, -10, 0, 10, mp.inf])
for e in ['0.2', '-0.2']:
    v = F_direct(mpf(e))
    print(f"F({e}) = {mp.nstr(v, 60)}")

# ---- Pade [N/D] for B(t) = sum b_n t^n, b_n = f_n/n! ----
NB = 61
b = [Fcoef[n]/factorial(n) for n in range(NB)]
def pade_coeffs(a, L, M):
    # solve for denominator coeffs q_1..q_M: sum_{j=0}^{M} q_j a_{L-M+k+j} = 0 ...
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            A[i, j] = a[L - M + i + j + 1] if L - M + i + j + 1 >= 0 else mpf(0)
        rhs[i] = -a[L + i + 1]
    q = (A**-1)*rhs
    den = [mpf(1)] + [q[i] for i in range(M)]
    num = [mpf(0)]*(L+1)
    for n in range(L+1):
        num[n] = sum(den[j]*a[n-j] for j in range(0, min(M, n)+1))
    return num, den
L, M = 30, 30
num, den = pade_coeffs(b, L, M)
def PB(t):
    nv = sum(num[k]*t**k for k in range(len(num)))
    dv = sum(den[k]*t**k for k in range(len(den)))
    return nv/dv
# pole check
import numpy as np
coefs = [complex(den[k]) for k in range(M, -1, -1)]
roots = np.roots(coefs)
im_poles = sorted([r for r in roots if abs(r.real) < 0.05], key=lambda r: abs(r.imag))[:6]
print("\nnearest-to-imaginary-axis poles of Pade[30/30] of B(t):")
for r in im_poles: print("   ", r, "   |t| =", abs(r))
offaxis = [r for r in roots if abs(r.real) >= 0.05]
print("poles off axis (count):", len(offaxis), " min |real| among them:", min(abs(r.real) for r in offaxis) if offaxis else None)
print("\nPade residual on series: max |num/den series mismatch| at order", NB, ":",
      mp.nstr(max(abs(sum(num[j]*b[n-j] for j in range(min(M,n)+1)) - (sum(den[j]*b[n-j] for j in range(min(M,n)+1)) if n<=L else 0)) for n in range(0, 5)), 3))
