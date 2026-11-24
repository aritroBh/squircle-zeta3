"""Numerical tests of Lemma B (branch), Lemma H (exact Hankel decomposition + uniform constants), Lemma R by hand.
F_t(z)=t1^(-2z/l) e^{(1+w)z} P_w(z)^(1/2), P_w(z)=1+e^{-2wz}+e^{-2z}; region III: w in [0,w0], w0=0.19, rho1=2."""
from mpmath import mp, mpf, mpc, exp, log, sqrt, pi, findroot, cos, sin, acos, beta, gamma, nstr, quad, re, im, arg, fabs
import numpy as np
mp.dps = 30
R = sqrt(pi**2+log(2)**2)/2; z0 = mpc(-log(2), pi)/2; th = arg(z0); rho1 = mpf(2); w0 = mpf('0.19')
P = lambda w, z: 1+exp(-2*w*z)+exp(-2*z)
dP = lambda w, z: -2*w*exp(-2*w*z)-2*exp(-2*z)
def zstar(w): return findroot(lambda z: P(w, z), z0*(1+w/2))
print("== Lemma R by hand: bound for min_{|z|=2}|2+e^{-2z}| in the band 1/2<r<7/2")
c1, c2 = -log(mpf(7)/2)/4, log(2)/4; s_min = min(sin(acos(c2)), sin(acos(c1)))
cpsi = cos(4*s_min)
m2 = 4-(2*cpsi)**2  # min_r 4+r^2+4r cos psi = 4-4cos^2psi
print(" sin(phi) >= %s, cos(psi) <= %s, |P_0|^2 >= %s, |P_0| >= %s ; e^{4 w0}-1 = %s" % (nstr(s_min, 8), nstr(cpsi, 8), nstr(m2, 8), nstr(sqrt(m2), 8), nstr(exp(4*w0)-1, 8)))
print("== Lemma B: branch of D_w^{1/2} along [0,z*], D_w(z)=P_w(z)/(1-z/z*)")
def track(w, zs, path):
    """continuous sqrt of D along path (list of points), start sqrt3"""
    prev = sqrt(3); vals = []
    for z in path:
        D = P(w, z)/(1-z/zs) if abs(z-zs) > mpf(10)**-20 else -zs*dP(w, zs)
        s = sqrt(D)
        if abs(s-prev) > abs(-s-prev): s = -s
        vals.append(s); prev = s
    return vals
for w in [mpf(0), mpf('0.05'), mpf('0.19')]:
    zs = zstar(w); path = [zs*mpf(j)/4000 for j in range(4001)]
    v = track(w, zs, path)
    # minimal Im along the way, and whether principal branch agrees
    pr = sqrt(-zs*dP(w, zs))
    imax = max(im(P(w, z)) for z in path[1:-1])
    print(" w=%s: z*=%s |z*|/R=%s  D^(1/2)(z*) tracked=%s principal=%s  max Im P_w(s z*) on (0,1)=%s" % (
        nstr(w, 3), nstr(zs, 12), nstr(abs(zs)/R, 10), nstr(v[-1], 12), nstr(pr, 12), nstr(imax, 4)))
psi0 = exp(z0)*sqrt(-z0*dP(0, z0))
print(" Psi_0 = e^{z0}(-z0 P_0'(z0))^{1/2} = %s ;  sqrt2*z0^(1/2) = %s" % (nstr(psi0, 15), nstr(sqrt(2)*sqrt(z0), 15)))
print("== Lemma H: exact decomposition [z^k]F = -(1/pi) 2Re[z*^-k J_k] + C_k  and leading-term error")
K = 60
def taylor_F(w, l, t1):
    # series of exp(a z) and sqrt(P_w)
    def expser(a): 
        out = [mpf(1)]
        for n in range(1, K+1): out.append(out[-1]*a/n)
        return out
    e1, e2 = expser(-2*w), expser(-2)
    Pser = [e1[n]+e2[n] for n in range(K+1)]; Pser[0] += 1
    q = [sqrt(Pser[0])]+[mpf(0)]*K
    for n in range(1, K+1): q[n] = (Pser[n]-sum(q[j]*q[n-j] for j in range(1, n)))/(2*q[0])
    pre = expser((1+w)-2*log(t1)/l)
    return [sum(pre[j]*q[n-j] for j in range(n+1)) for n in range(K+1)]
xg, wg = np.polynomial.legendre.leggauss(160)
def gl(f, a, b):
    return sum(mpf(wi)*(b-a)/2*f((b-a)/2*mpf(xi)+(b+a)/2) for xi, wi in zip(xg, wg))
def D2sqrt(w, zs, z):
    # D2 = P/((1-z/zs)(1-z/conj zs)) zero-free on |z|<=2; its sqrt: track radially from 0
    n = 200; prev = sqrt(3)
    for j in range(1, n+1):
        u = z*mpf(j)/n
        D = P(w, u)/((1-u/zs)*(1-u/zs.conjugate()))
        s = sqrt(D)
        if abs(s-prev) > abs(-s-prev): s = -s
        prev = s
    return prev
for (w, l) in [(mpf(0), mpf(40)), (mpf('0.05'), mpf(3)), (mpf('0.19'), mpf(1)), (mpf('0.1'), mpf(10))]:
    lam = w*l; t1 = 1/(1+exp(-lam)+exp(-l)); zs = zstar(w); v1 = rho1/abs(zs)-1
    Fk = taylor_F(w, l, t1)
    pref = lambda z: exp(-2*z*log(t1)/l)*exp((1+w)*z)
    # check: Re D2 > 0 on disc? sample
    def Psi(z): return pref(z)*sqrt(1-z/zs.conjugate())*D2sqrt(w, zs, z)
    def Fcirc(phi):
        z = rho1*exp(mpc(0, phi))
        return pref(z)*sqrt(1-z/zs)*sqrt(1-z/zs.conjugate())*D2sqrt(w, zs, z)
    a1 = arg(zs)
    psivals = {}
    for k in [5, 10, 20, 40, 60]:
        J = gl(lambda u: 2*u**2*Psi(zs*(1+u**2))*(1+u**2)**(-k-1), mpf(0), sqrt(v1))
        cut = -2/pi*re(zs**(-k)*J)
        circ = sum(gl(lambda ph: Fcirc(ph)*(rho1*exp(mpc(0, ph)))**(-k), a, b) for a, b in [(-pi, -a1), (-a1, a1), (a1, pi)])/(2*pi)
        lead = -2/pi*re(zs**(-k)*Psi(zs))*beta(mpf(3)/2, k-mpf(1)/2)
        print(" w=%s l=%s k=%2d: [z^k]F=%s  identity resid=%s  (lead-exact)/scale=%s  circ/scale=%s" % (
            nstr(w, 3), nstr(l, 3), k, nstr(Fk[k], 12), nstr(abs(cut+re(circ)-Fk[k])/abs(zs)**(-k), 3),
            nstr((lead-Fk[k])/(abs(zs)**(-k)*k**-1.5), 4), nstr(abs(circ)/(abs(zs)**(-k)*k**-1.5), 3)))
