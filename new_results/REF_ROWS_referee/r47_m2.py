# R47: M_2(alpha) from the corner integral (Cauchy integral in u), vs Thm M2alpha (HPL form, own HPL series)
# and vs the classical form B2' typed from research/10_frontier/G1_zeta3_any_angle/result.md.
from mpmath import mp, mpf, mpc, sqrt, log, cos, sin, pi, quad, zeta, polylog, nstr, exp, atanh, fsum
import sys
mp.dps = 50
NPTS = 64; RAD = mpf('0.3')
def Mu(u, ca):
    f = lambda g: ((1+g)**(-1-u)*sqrt(1+g**(2-2*u)-2*ca*g**(1-u))-1)*g**(u-1)
    return quad(f, [0, mpf(1)/10**6, mpf(1)/1000, mpf(1)/10, mpf(1)/2, 1])
def Mk(alpha, K=3):
    ca = cos(alpha)
    vals = []
    for j in range(NPTS//2+1):
        u = RAD*exp(2j*pi*j/NPTS)
        vals.append(Mu(u, ca))
    full = vals + [v.conjugate() for v in vals[1:NPTS//2][::-1]]
    out = []
    for k in range(K+1):
        s = fsum(full[j]*exp(-2j*pi*j*k/NPTS) for j in range(NPTS))/NPTS/RAD**k
        out.append(s.real)
    return out
# HPL series for words in {0,1,-1} with last letter nonzero
def hpl(word, x, nmax):
    # a[n] coefficients of x^n, n=0..nmax
    w = list(word)
    last = w.pop()
    a = [mpf(0)]*(nmax+1)
    for n in range(1, nmax+1):
        a[n] = mpf(1)/n if last == 1 else mpf((-1)**(n-1))/n
    for letter in reversed(w):
        b = [mpf(0)]*(nmax+1)
        if letter == 0:
            for n in range(1, nmax+1): b[n] = a[n]/n
        elif letter == 1:
            cum = mpf(0)
            for n in range(1, nmax+1):
                cum += a[n-1]; b[n] = cum/n
        else:
            cum = mpf(0)  # coefficient of t^{n-1} in a(t)/(1+t) = sum_{j<=n-1} a_j (-1)^{n-1-j}
            for n in range(1, nmax+1):
                cum = -cum + a[n-1]; b[n] = cum/n
        a = b
    return fsum(a[n]*x**n for n in range(1, nmax+1))
def Y(x, nmax):
    H = lambda w: hpl(w, x, nmax)
    return (1-x)*(H((0,0,1))+H((0,1,-1))-H((1,-1,1))) + (1-x**2)/(2*x)*H((1,0,1)) - (1-x)**2/(2*x)*H((1,0,-1))
def Phi(c): return mpf(1)/2 - (1-c**2)/(2*c)*atanh(c)
def M2_hpl(c):
    nmax = int(mp.dps*2.4/(-log(c, 10)))+50
    return zeta(2)*Phi(c) - polylog(2, c**2)/2 + Y(c, nmax) + Y(-c, nmax)
def M2_classical(c):
    lp = log(1+c); lm = log(1-c); l0 = log(c); lam = log(2); q = (1-c**2)/(2*c)
    L2 = lambda z: polylog(2, z); L3 = lambda z: polylog(3, z)
    return (zeta(2)*Phi(c) - L2(c) - L2(-c)
      + (1+c)*L3(c) + (1-c)*L3(-c) + (1-c)/c*L3(1-c) - (1+c)/c*L3(1/(1+c)) - 2*c*L3((1+c)/2)
      + (1-c)*(L3((1-c)/(1+c)) + L3(-(1-c)/(1+c))) - 2*c*L3(2*c/(1+c))
      + L2((1+c)/2)*((1+c)*lp - (1-c)*lm) + q*(L2(c)-L2(-c))*(lp+lm)
      - (1-c)*lm**2*lp + (1-c)*lm*lp**2 + (4*c**2-c+1)/(6*c)*lp**3 + q*l0*lm**2 - (1-c)*l0*lm*lp - (1+c**2)/(2*c)*l0*lp**2
      + lam*((1-c)*lm**2 - (1-c)*lm*lp - c*lp**2) + lam**2*(1-c)/2*(lp-lm) - pi**2*((1-c)/(6*c)*lm + (1+c)/(6*c)*lp)
      + (mpf(7)/4 + 2*c)*zeta(3) - c/6*pi**2*lam + c/3*lam**3)
if __name__ == "__main__":
    for alpha in [pi/3, pi/2, 2*pi/3, mpf('2.6'), mpf('1.0')]:
        c = cos(alpha/2)
        mk = Mk(alpha)
        print("alpha", nstr(alpha, 6), " c", nstr(c, 8))
        M0 = -(1+c)*log(1+c)-(1-c)*log(1-c)
        print("  M0 check", nstr(abs(mk[0]-M0), 3))
        h = M2_hpl(c); cl = M2_classical(c)
        print("  M2 integral", nstr(mk[2], 35))
        print("  |integral - HPL form|", nstr(abs(mk[2]-h), 3), "  |integral - classical B2'|", nstr(abs(mk[2]-cl), 3), "  |HPL - classical|", nstr(abs(h-cl), 3))
    sys.stdout.flush()
