"""
angle.py -- Problem 3 at an arbitrary corner angle theta (combining Ahmed's Thm C with general f):
  gamma_n(x) = x u + f(x)^n v,  |u|=|v|=1,  u.v = c = cos(theta).
  D_n^theta = 2 - length = int_0^1 h_c(y_n'(x)) dx,  h_c(s) = 1 + s - sqrt(1 + 2cs + s^2).
Layer cake with density h_c'(s) = 1 - (s+c)/sqrt(1+2cs+s^2) (mass 1-c).
Claim (Corollary 5):  D_n^theta = [ (1-c) log(n f'(1)) + C(theta) ]/(f'(1)(n-1)) + O(log^2 n/n^2),
  C(theta) = 2 sin^2(theta/2) + 4 cos^2(theta/2) log cos(theta/2)   (Ahmed's constant),
and the full expansion has the same structure as for theta = pi/2, with the moments of log(1/S)
replaced by rho_m(c) = int_0^inf h_c'(s) log^m(1/s) ds (computed here by quadrature, 60 digits).
"""
import mpmath as mp
import sympy as sp
from fnlib import TESTS, make_callables, jets, solve_T

def hc(s, c):
    return 2*(1 - c)*s/(1 + s + mp.sqrt(1 + 2*c*s + s*s))

def hcp(s, c):
    r = mp.sqrt(1 + 2*c*s + s*s)
    return (1 - c*c)/(r*(r + s + c))       # = 1 - (s+c)/r, cancellation-free

def rho(c, M):
    out = []
    for m in range(M + 1):
        out.append(mp.quad(lambda u: hcp(mp.exp(u), c)*mp.exp(u)*(-u)**m, [-mp.inf, -20, -5, 0, 5, 20, mp.inf]))
    return out

def tail_series(c, kmax):
    """h_c'(s) = sum_{k>=2} d_k s^{-k} for large s"""
    t = sp.Symbol('t', positive=True)
    cc = sp.nsimplify(c) if isinstance(c, (int,)) else sp.Rational(str(c)) if not isinstance(c, sp.Basic) else c
    expr = 1 - (1/t + cc)/sp.sqrt(1 + 2*cc/t + 1/t**2)
    ser = sp.series(expr, t, 0, kmax + 1).removeO()
    return {k: mp.mpf(sp.N(ser.coeff(t, k), mp.mp.dps + 10)) for k in range(2, kmax + 1)}

def expansion_polys_c(A, B, K, c, rh, dk):
    a = A[1]
    T = solve_T(A, B, K)
    la = mp.log(a)
    P = [[mp.mpf(0)]*(k + 1) for k in range(K + 1)]
    for (i, j), t in T.items():
        k = i + j
        for m in range(i + 1):
            P[k][i - m] += t*mp.binomial(i, m)*rh[m]
    for (i, j), t in T.items():
        k0 = i + j
        for kk, dcoef in dk.items():
            p = kk - 1                     # (n a)^{1-kk} = a^{1-kk} mu^{p} (1+mu)^{-p}
            if k0 + p > K:
                continue
            sm = mp.mpf(0)
            for l in range(i + 1):
                sm += mp.binomial(i, l)*la**(i - l)*mp.factorial(l)/mp.mpf(kk - 1)**(l + 1)
            coef = dcoef*(-1)**i*sm*a**(1 - kk)
            for r in range(0, K - k0 - p + 1):
                P[k0 + p + r][0] -= t*coef*mp.binomial(-p, r)
    return P

def defect_c(f, fp, n, a, c, dps_target):
    n = mp.mpf(n); nm1 = n - 1
    vmax = min((dps_target + 15)*mp.log(10)/a + 2*mp.log(n*a)/a + 50, n)
    def integrand(v):
        xx = 1 - v/n
        fx = f(xx)
        if fx <= 0:
            return mp.mpf(0)
        return hc(n*mp.exp(nm1*mp.log(fx))*fp(xx), c)
    cen = mp.log(n*a)/a
    pts = sorted(set([mp.mpf(0), vmax] + [cen + mp.mpf(k)/a for k in [-8, -4, -2, -1, 0, 1, 2, 4, 8, 16, 32, 64, 128, 256]
                                          if 0 < cen + mp.mpf(k)/a < vmax]))
    return mp.quad(integrand, pts)/n

if __name__ == '__main__':
    K = 9
    for (theta_name, cval) in [('pi/3', sp.Rational(1, 2)), ('2pi/3', sp.Rational(-1, 2))]:
        mp.mp.dps = 110
        c = mp.mpf(sp.N(cval, 120))
        theta = mp.acos(c)
        Cth = 2*mp.sin(theta/2)**2 + 4*mp.cos(theta/2)**2*mp.log(mp.cos(theta/2))
        rh = rho(c, K)
        dk = tail_series(cval, K + 2)
        print('=' * 100)
        print('theta = %s:  rho_0 - (1-c) = %s ;  rho_1 - C(theta) = %s   (Ahmed\'s C(theta) = %s)' % (
            theta_name, mp.nstr(rh[0] - (1 - c), 3), mp.nstr(rh[1] - Cth, 3), mp.nstr(Cth, 30)))
        for name in ['expm1', 'mobius']:
            expr = TESTS[name]
            f, fp = make_callables(expr)
            A, B = jets(expr, K + 2)
            a = A[1]
            P = expansion_polys_c(A, B, K, c, rh, dk)
            Cf = (1 - c)*mp.log(a) + Cth
            print('  f = %-7s exact constant (1-c) log f\'(1) + C(theta) = %s' % (name, mp.nstr(Cf, 30)))
            print('      P_1 from engine: [%s, %s]   expected [C/a, (1-c)/a] = [%s, %s]' % (
                mp.nstr(P[1][0], 20), mp.nstr(P[1][1], 20), mp.nstr(Cf/a, 20), mp.nstr((1 - c)/a, 20)))
            Ds = {}
            for e in range(6, 21):
                n = 10**e
                Ds[n] = defect_c(f, fp, n, a, c, 90)
            # single-n extraction with engine K
            for e in [6, 10, 15, 20]:
                n = 10**e; mu = 1/mp.mpf(n - 1); L = mp.log(n)
                hi = sum(mu**k*mp.polyval(P[k][::-1], L) for k in range(2, K + 1))
                Cest = (a/mu)*(Ds[n] - hi) - (1 - c)*L
                print('      n=1e%-2d  C_est - C = %s' % (e, mp.nstr(Cest - Cf, 3)))
            # theory-free fit, unknown P_1..P_5 (20 unknowns), n = 1e6..1e20 (15 points) + fixed tail P_6..P_9
            cols = [(k, i) for k in range(1, 5) for i in range(k + 1)]   # 14 unknowns
            nsel = sorted(Ds)[-len(cols):]
            Mx = mp.matrix(len(cols), len(cols)); rhs = mp.matrix(len(cols), 1)
            for r, n in enumerate(nsel):
                mu = 1/mp.mpf(n - 1); L = mp.log(n)
                for j, (k, i) in enumerate(cols):
                    Mx[r, j] = mu**k*L**i
                rhs[r] = Ds[n] - sum(mu**k*mp.polyval(P[k][::-1], L) for k in range(5, K + 1))
            sol = mp.lu_solve(Mx, rhs)
            print('      Richardson (P_1..P_4 unknown, 14 coeffs; P_5..P_9 from engine): C_fit - C = %s' % mp.nstr(sol[0]*a - Cf, 3))
            print()
