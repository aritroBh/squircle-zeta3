
import json, os
from mpmath import mp, mpf, mpc, pi, factorial, quad, sqrt, cosh, exp, log, sin, cos, nstr

DIR = os.path.dirname(os.path.abspath(__file__))

def load_f(nmax=64):
    with open(os.path.join(DIR, 'step2_num200.json')) as fh:
        raw = json.load(fh)['f200']
    return [mpf(v) for v in raw[:nmax+1]]

def borel_coeffs(f):
    return [f[n]/factorial(n) for n in range(len(f))]

def pade_coeffs(a, L, M):
    # solve the [L/M] Pade of series sum a_n t^n (standard Hankel system)
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(M):
        for j in range(M):
            k = L + i - j
            A[i, j] = a[k] if k >= 0 else mpf(0)
        rhs[i] = -a[L + 1 + i]
    q = (A**-1)*rhs
    den = [mpf(1)] + [q[j] for j in range(M)]
    num = [sum(den[j]*a[n-j] for j in range(0, min(M, n)+1)) for n in range(L+1)]
    return num, den

def pade_eval(num, den, t):
    nv = num[0]
    for k in range(1, len(num)): nv += num[k]*t**k
    dv = den[0]
    for k in range(1, len(den)): dv += den[k]*t**k
    return nv/dv

def F_direct(eta, dps_quad=60):
    # direct quadrature of F(eta) = int_0^1 (t(1-t))^{eta/2} sqrt(cosh(eta*l)) dt, l=log(t/(1-t))
    old = mp.dps; mp.dps = dps_quad
    try:
        f = lambda t: (t*(1-t))**(eta/2) * sqrt(cosh(eta*log(t/(1-t))))
        v = quad(f, [0, mpf('0.25'), mpf('0.5'), mpf('0.75'), 1])
    finally:
        mp.dps = old
    return v

# ---- conformal map: t-plane cut along +-i[pi/2, i inf) -> unit w-disc ----
def w_of_t(t):
    # u = 2t/(i*pi) = -2i t/pi ; u = 2w/(1+w^2) maps disc to C \ ((-inf,-1] u [1,inf))
    u = -2j*t/pi if isinstance(t, complex) else mpc(0, -2)*t/pi
    u = mpc(u)
    return (1 - sqrt(1 - u*u))/u

def t_of_w(w):
    return mpc(0, pi)*w/(1+w*w)
