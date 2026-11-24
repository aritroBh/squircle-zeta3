# P2N4R referee: own numerical check of Lemma g3:prop:n4 from the exact formula of thm:P2 (mpmath, 110 digits).
from mpmath import mp, mpf, mpc, pi, sin, cos, log, exp, hyp2f1, polylog, zeta, quad, nstr, fabs, expj, sqrt, legendre
import sys, time
mp.dps = 110
t0 = time.time()
def lam_m_list(c, M):                      # Lambda_m = (P_{m-1}-P_{m+1})/(2m+1), Legendre by three-term recurrence
    P = [mpf(1), c]
    for j in range(1, M+2): P.append(((2*j+1)*c*P[j] - j*P[j-1])/(j+1))
    return [None] + [(P[m-1]-P[m+1])/(2*m+1) for m in range(1, M+1)]
def g(e, z):                               # R_theta(eps)/(2z) = pi e/sin(pi e) 2F1(-e,1+e;2;z)
    return pi*e/sin(pi*e)*hyp2f1(-e, 1+e, 2, z)
def Lexact(th, n):                          # Theorem 15 (thm:P2)
    c = cos(th); z = sin(th/2)**2; n = mpf(n); e = 1/(n-1)
    M = int((mp.dps+10)*log(10)/log(n)) + 3
    Lm = lam_m_list(c, M)
    s = 1 + c + 2*z*g(e, z)*n**(-e)
    for m in range(1, M+1): s += (-1)**m*e/(m-e)*Lm[m]*n**(-m)
    return s
def closed(z):
    w = (1-z)/z; l = log(1-z); u = w*l
    S12 = lambda t: zeta(3) - polylog(3, 1-t) + log(1-t)*polylog(2, 1-t) + log(1-t)**2*log(t)/2
    S22 = quad(lambda t: S12(t)/t, [0, z/4, z/2, z])
    return (-1-u, zeta(2)+mpf(1)/2-polylog(2, z)-u**2/2, 2*w*S12(z)-(1+u**3)/3,
            zeta(4)/2 - 2*S22 + 2*w*w*l*S12(z) + (1-u**4)/4)
def cauchy(z, J=5, N=256, r=mpf(1)/5):
    vals = [mp.log(g(r*expj(2*pi*(k+mpf(1)/2)/N), z)) for k in range(N)]
    return [(sum(vals[k]*expj(-2*pi*(k+mpf(1)/2)*j/N) for k in range(N))/N/r**j).real for j in range(1, J+1)]
print("== (0) Lexact vs direct quadrature of the arc length (40 digits)"); sys.stdout.flush()
mp.dps = 45
for th in (mpf('0.4'), mpf('2.9')):
    for n in (9, mpf('13.5')):
        c = cos(th); nn = mpf(n)
        f = lambda t: sqrt(1 + 2*c*nn*t**(nn-1) + nn*nn*t**(2*nn-2))
        val = quad(f, [0, mpf(1)/2, mpf(3)/4, mpf(7)/8, mpf(15)/16, 1])
        print("  th=%s n=%s  quad - exact = %s" % (nstr(th, 3), nstr(nn, 4), nstr(val - Lexact(th, nn), 3)))
mp.dps = 110
angles = [mpf('1e-3'), mpf('0.3'), mpf(1), pi/2, mpf('2.5'), pi-mpf('0.01'), pi-mpf('1e-3')]
ns = [9, mpf('9.5'), 10, 20, 50, 100, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
sup1 = sup2 = sup3 = 0
print("== (1) residual R_n = E_n - (d-2z a2 L)/n^3 - Q4/n^4.")
print("   col1 = R n^5/(z L^3)   [lemma (b): |col1| <= 1.2e5; limit -a2/3]")
print("   col2 = (R + z a2 L^3/(3n^5)) n^5/(z L^2)   [lemma (a): |col2| <= 2.5e5]")
print("   col3 = (R n^5 - Q5(L)) n/(z L^4), Q5 = full formal x^5 coefficient (needs a5)  [bounded <=> Q5 is the next term]")
for th in angles:
    z = sin(th/2)**2
    a1, a2, a3, a4 = closed(z)
    cf = cauchy(z)
    diff = max(fabs(p-q) for p, q in zip((a1, a2, a3, a4), cf[:4]))
    a5 = cf[4]
    c2 = 2*z*(a2-1+z); d = 2*z*(a2*(2+a1)+a3) - z*(3-z-2*z*z)
    p4 = -z*(1-z)*(20*z*z-2*z+43)/6
    print(" th=%s z=%s  a1..a4 closed vs Cauchy: max diff %s ; a2=%s a3=%s a4=%s a5=%s ; -a2/3=%s  (%.0fs)" % (
        nstr(th, 6), nstr(z, 6), nstr(diff, 3), nstr(a2, 8), nstr(a3, 8), nstr(a4, 8), nstr(a5, 8), nstr(-a2/3, 8), time.time()-t0))
    for n in ns:
        n = mpf(n); e = 1/(n-1); L = log(n); lam = L - a1
        E = Lexact(th, n) - (1 + cos(th) + 2*z*exp(a1*e)*n**(-e) + c2/n**2)
        Q4 = z*(a2*((lam-3)**2-3) - 2*a3*(lam-3) + a2**2 + 2*a4) + p4
        Q5 = (-a2*z/3*L**3 + z*(a1*a2+4*a2+a3)*L**2 - z*(a1**2*a2+8*a1*a2+2*a1*a3+a2**2+12*a2+8*a3+2*a4)*L
              + (a1**3*a2*z/3 + 4*a1**2*a2*z + a1**2*a3*z + a1*a2**2*z + 12*a1*a2*z + 8*a1*a3*z + 2*a1*a4*z + 4*a2**2*z
                 + 2*a2*a3*z + 8*a2*z + 12*a3*z + 8*a4*z + 2*a5*z + 7*z**5 - mpf(235)/18*z**4 + mpf(191)/18*z**3
                 + mpf(115)/12*z**2 - mpf(509)/36*z))
        R = E - (d - 2*z*a2*L)/n**3 - Q4/n**4
        c1 = R*n**5/(z*L**3); c2_ = (R + z*a2*L**3/(3*n**5))*n**5/(z*L**2); c3 = (R*n**5 - Q5)*n/(z*L**4)
        sup1 = max(sup1, fabs(c1)); sup2 = max(sup2, fabs(c2_)); sup3 = max(sup3, fabs(c3))
        print("   n=%-8s R=%-11s col1=%-11s col2=%-11s col3=%s" % (nstr(n, 3), nstr(R, 4), nstr(c1, 6), nstr(c2_, 6), nstr(c3, 5)))
    sys.stdout.flush()
print("sup|col1| = %s (bound 1.2e5) ; sup|col2| = %s (bound 2.5e5) ; sup|col3| = %s" % (nstr(sup1, 5), nstr(sup2, 5), nstr(sup3, 5)))
mp.dps = 30
print("== (2) sup |H| on |eps| = 1/4 (the 3.1 of the Cauchy bounds), z grid, 96 points")
best = 0
for z in [mpf(k)/40 for k in range(1, 40)] + [mpf('0.001'), mpf('0.999')]:
    _, a2, _, _ = closed(z); a1 = -1-(1-z)/z*log(1-z)
    for k in range(96):
        e = expj(2*pi*k/96)/4
        best = max(best, fabs(g(e, z)*exp(-a1*e)))
print("  max |H| =", nstr(best, 6), "(proof uses 3.1)")
print("done in %.0fs" % (time.time()-t0))
