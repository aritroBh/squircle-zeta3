"""step_01.py -- derive the w-integrand for the corner family and check numerically.

CLAIM (derived in result.md): with t = g/(1+g), the corner function is
  M(u;alpha) = int_0^1 Psi_u(t) dt,  Psi_u = (sqrt(t^{2u-2} + (1-t)^{2u-2} - 2 ca (t(1-t))^{u-1}) - t^{u-1} - (1-t)^{u-1})/2,
ca = cos(alpha). With c = cos(alpha/2), r = (1-c)/(1+c) and t = (w-r)(w+1)/(2(1-r)w):
  M(u;alpha) = int_r^1 Psi_u(t(w)) (dt/dw) dw,   dt/dw = (w^2+r)/(2(1-r) w^2),
and the master sqrt rationalises: 1 - 4 c^2 t(1-t) = (w^2+r)^2/((1+r)^2 w^2).

Checks vs corner_series.py (M_exact and corner_coeffs) at two angles.
"""
import sys
sys.path.insert(0, '.')
import mpmath as mp
mp.mp.dps = 50
from corner_series import corner_coeffs, M_exact

def psi(t, u, ca):
    return (mp.sqrt(t**(2*u-2) + (1-t)**(2*u-2) - 2*ca*(t*(1-t))**(u-1)) - t**(u-1) - (1-t)**(u-1))/2

def M_w(u, alpha):
    c = mp.cos(alpha/2); ca = mp.cos(alpha); r = (1-c)/(1+c)
    def t_of(w): return (w-r)*(w+1)/(2*(1-r)*w)
    def f(w):
        return psi(t_of(w), u, ca) * (w*w+r)/(2*(1-r)*w*w)
    return mp.quad(f, [r, mp.sqrt(r), 1])

def M_t(u, alpha):
    ca = mp.cos(alpha)
    return mp.quad(lambda t: psi(t, u, ca), [0, mp.mpf('0.5'), 1])

for alpha in (mp.pi*2/3, mp.pi/3):
    for u in (mp.mpf('0.3'), mp.mpf('-0.25')):
        a = M_exact(u, alpha); b = M_t(u, alpha); cc = M_w(u, alpha)
        print(f"alpha={mp.nstr(alpha,6)} u={u}: exact={mp.nstr(a,20)}")
        print(f"   t-int={mp.nstr(b,20)} diff={mp.nstr(abs(a-b),3)} ; w-int={mp.nstr(cc,20)} diff={mp.nstr(abs(a-cc),3)}")

# Taylor coefficients: M_k(alpha) = [u^k] int Psi.  Use the series of the w-integrand.
def corner_coeffs_w(alpha, N):
    c = mp.cos(alpha/2); ca = mp.cos(alpha); r = (1-c)/(1+c)
    def t_of(w): return (w-r)*(w+1)/(2*(1-r)*w)
    dt = lambda w: (w*w+r)/(2*(1-r)*w*w)
    # series in u at each w, exact formal expansion like corner_series
    def series_at(w, N):
        t = t_of(w)
        L1 = mp.log(t); L2 = mp.log1p(-t)
        fac = [mp.factorial(k) for k in range(N+1)]
        E1 = [L1**k/fac[k] for k in range(N+1)]      # t^u
        E2 = [L2**k/fac[k] for k in range(N+1)]      # (1-t)^u
        # sqrt(t^{-2}t^{2u} + (1-t)^{-2}(1-t)^{2u} - 2ca t^{-1}(1-t)^{-1} (t(1-t))^u)
        Q = [t**-2 + (1-t)**-2 - 2*ca*t**-1*(1-t)**-1]
        for n in range(1, N+1):
            Q.append(t**-2*(2*L1)**n/fac[n] + (1-t)**-2*(2*L2)**n/fac[n] - 2*ca*t**-1*(1-t)**-1*(L1+L2)**n/fac[n])
        # sqrt series
        B = [mp.sqrt(Q[0])]
        for n in range(1, N+1):
            s = Q[n] - mp.fsum(B[k]*B[n-k] for k in range(1, n))
            B.append(s/(2*B[0]))
        out = []
        for n in range(N+1):
            val = B[n] - (t**-1*E1[n] + (1-t)**-1*E2[n])
            out.append(val/2 * dt(w))
        return out
    cache = {}
    def F(w, k):
        if w not in cache: cache[w] = series_at(w, N)
        return cache[w][k]
    return [mp.quad(lambda w: F(w, k), [r, mp.sqrt(r), 1]) for k in range(N+1)]

for alpha in (mp.pi*2/3, mp.pi/3):
    Mk = corner_coeffs(alpha, 4)
    Mkw = corner_coeffs_w(alpha, 4)
    print("alpha =", mp.nstr(alpha, 8))
    for k in range(5):
        print(f"  M_{k}: corner={mp.nstr(Mk[k],25)}  w={mp.nstr(Mkw[k],25)}  diff={mp.nstr(abs(Mk[k]-Mkw[k]),3)}")
print("STEP01 DONE")
