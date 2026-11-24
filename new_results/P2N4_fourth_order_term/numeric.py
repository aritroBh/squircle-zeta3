# P2N4: numerical check of the n^-4 lemma against the exact formula of thm:P2 (mpmath, 90 digits).
from mpmath import (mp, mpf, hyp2f1, sin, cos, pi, log, exp, polylog, zeta, quad, nstr, fabs, expj, sqrt)
mp.dps = 90
def S12(z):
    l = log(1-z); return zeta(3) - polylog(3, 1-z) + l*polylog(2, 1-z) + l*l*log(z)/2
def S22(z): return quad(lambda t: S12(t)/t, [0, z/2, z])
def avals(z):
    w = (1-z)/z; l = log(1-z); u = w*l
    a1 = -1-u; a2 = zeta(2)+mpf(1)/2-polylog(2, z)-u**2/2; a3 = 2*w*S12(z)-(1+u**3)/3
    a4 = zeta(4)/2 - 2*S22(z) + 2*w*w*l*S12(z) + (1-u**4)/4
    return a1, a2, a3, a4
def g(e, z): return pi*e/sin(pi*e)*hyp2f1(-e, 1+e, 2, z)
def Lex(th, n):   # exact formula of thm:P2
    c = cos(th); z = sin(th/2)**2; n = mpf(n); e = 1/(n-1)
    s = 1 + c + 2*z*g(e, z)*n**(-e)
    P = [mpf(1), c]; m = 1
    while True:
        while len(P) < m+2:
            j = len(P)-1; P.append(((2*j+1)*c*P[j]-j*P[j-1])/(j+1))
        s += (-1)**m*e/(m-e)*(P[m-1]-P[m+1])/(2*m+1)*n**(-m)
        if 4*n**(-m) < mpf(10)**(-mp.dps-5): break
        m += 1
    return s
# (0) a_1..a_4 closed forms vs Taylor coefficients of log g by a Cauchy integral (radius 1/5, 160 nodes)
print("== a1..a4 closed forms vs Cauchy-integral coefficients of log g (r=1/5, N=160)")
worst = 0
for th in (mpf('0.3'), pi/2, mpf('2.8')):
    z = sin(th/2)**2; N = 160; r = mpf(1)/5
    vals = [mp.log(g(r*expj(2*pi*(k+mpf(1)/2)/N), z)) for k in range(N)]
    cf = [(sum(vals[k]*expj(-2*pi*(k+mpf(1)/2)*j/N) for k in range(N))/N/r**j).real for j in range(1, 5)]
    worst = max([worst] + [fabs(a-b) for a, b in zip(avals(z), cf)])
print("  worst difference:", nstr(worst, 3), "(aliasing error of the rule is about (1/5/(1/4))^160 ~ 1e-15.5 relative to |a_{j+160}|r^160; here tiny)")
# (1) independent check of the exact formula by direct quadrature of the arc length
print("== direct quadrature of L_theta(n) vs exact formula (50 digits)")
mp.dps = 50
for th in (mpf('0.7'), mpf('2.4')):
    for n in (50, 400):
        c = cos(th); nn = mpf(n)
        f = lambda t: sqrt(1+2*c*nn*t**(nn-1)+nn*nn*t**(2*nn-2))
        pts = sorted(set([mpf(0)]+[1-mpf(j)/nn for j in (160, 80, 40, 20, 10, 5, 2, 1) if j < nn]+[mpf(1)]))
        print("  th=%s n=%d: direct - exact = %s" % (nstr(th, 3), n, nstr(quad(f, pts, maxdegree=12)-Lex(th, n), 3)))
mp.dps = 90
# (2) the residual R_n = E_n - (d - 2 z a2 L)/n^3 - Q4/n^4
print("== residual R_n after Q4.  col1 = R n^5/(z log^3 n) (lemma: |.|<=C4); col2 = (R + z a2 log^3n/(3n^5)) n^5/(z log^2 n) (lemma: |.|<=C')")
print("   limit of col1 as n->oo is -a2/3")
ns = [9, 10, 20, 50, 100, 10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
s1 = s2 = 0
for th in (mpf('0.01'), mpf('0.3'), mpf(1), pi/2, mpf(2), mpf('2.8'), mpf('3.13')):
    z = sin(th/2)**2; a1, a2, a3, a4 = avals(z)
    c2 = 2*z*(a2-1+z); d = 2*z*(a2*(2+a1)+a3) - z*(3-z-2*z*z)
    p4 = -z*(1-z)*(20*z*z-2*z+43)/6
    print(" th=%s z=%s  limit col1 = -a2/3 = %s ; formal limit col2 = a1 a2 + 4 a2 + a3 = %s" % (nstr(th, 4), nstr(z, 4), nstr(-a2/3, 6), nstr(a1*a2+4*a2+a3, 6)))
    for n in ns:
        n = mpf(n); e = 1/(n-1); L = log(n); lam = L - a1
        E = Lex(th, n) - (1 + cos(th) + 2*z*exp(a1*e)*n**(-e) + c2/n**2)
        Q4 = z*(a2*((lam-3)**2-3) - 2*a3*(lam-3) + a2**2 + 2*a4) + p4
        R = E - (d-2*z*a2*L)/n**3 - Q4/n**4
        r1 = R*n**5/(z*L**3); r2 = (R + z*a2*L**3/(3*n**5))*n**5/(z*L**2)
        s1 = max(s1, fabs(r1)); s2 = max(s2, fabs(r2))
        print("   n=%-9s R=%-12s col1=%-12s col2=%s" % (nstr(n, 3), nstr(R, 4), nstr(r1, 6), nstr(r2, 6)))
print("sup |col1| =", nstr(s1, 5), " (proof: C4 = 1.2e5);  sup |col2| =", nstr(s2, 5), " (proof: C' = 2.5e5)")
