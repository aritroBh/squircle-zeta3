# R36: M_1(alpha) from the corner integral vs Prop. M1short and the printed 2pi/3 form.
from mpmath import mp, mpf, sqrt, log, cos, pi, polylog, quad, zeta, nstr, diff
mp.dps = 60
def M1_def(alpha):
    ca = cos(alpha)
    def f(g):
        Q = sqrt(1+g*g-2*ca*g)
        # d/du at u=0 of [(1+g)^{-1-u} sqrt(1+g^{2-2u}-2ca g^{1-u}) - 1] g^{u-1}
        t1 = (-log(1+g)*Q + (g-ca)*(-g*log(g))/Q)/(1+g)/g
        t2 = (Q/(1+g)-1)*log(g)/g
        return t1+t2
    return quad(f, [0, mpf(1)/1000, mpf(1)/10, mpf(1)/2, 1])
def M_u(u, alpha):  # full corner integral, for a finite-difference check
    ca = cos(alpha)
    return quad(lambda g: ((1+g)**(-1-u)*sqrt(1+g**(2-2*u)-2*ca*g**(1-u))-1)*g**(u-1), [0, mpf(1)/1000, mpf(1)/10, 1])
def M0(c): return -(1+c)*log(1+c)-(1-c)*log(1-c)
def M1_prop(c):
    lam = log(2); L2 = lambda x: polylog(2, x)
    return zeta(2)-lam**2-lam*M0(c)-(1-c)*(L2(c)+L2((1-c)/2))-(1+c)*(L2(-c)+L2((1+c)/2))
print("Prop M1short vs integral definition")
for a in [mpf('0.3'), pi/3, pi/2, 2*pi/3, mpf('2.6'), mpf('3.0')]:
    c = cos(a/2)
    d = M1_def(a); p = M1_prop(c)
    print(nstr(a,6), nstr(d,25), "diff", nstr(abs(d-p),3))
# M_0 check against u=0 integral
a = 2*pi/3
print("M0(2pi/3) integral - formula:", nstr(abs(M_u(0,a)-M0(cos(a/2))),3))
# finite-difference check of the derivative formula itself
mp.dps = 40
h = mpf('1e-8')
fd = (M_u(h,a)-M_u(-h,a))/(2*h)
print("finite diff vs M1_def at 2pi/3:", nstr(abs(fd-M1_def(a)),3))
mp.dps = 120
lam = log(2); l3 = log(3)
printed = (pi**2+48*lam**2-6*l3**2-24*lam*l3-12*polylog(2, mpf(1)/3))/24
print("M1(2pi/3) printed =", nstr(printed, 40))
print("printed - Prop(c=1/2):", nstr(abs(printed-M1_prop(mpf(1)/2)),3))
mp.dps=80
print("printed - integral def (80 dps):", nstr(abs(printed-M1_def(2*pi/3)),3))
# square corner: Prop at c=1/sqrt2 vs Theorem M1
mp.dps=60
s2=sqrt(2); L=log(1+s2)
thm = (19*s2/48-mpf(1)/24)*pi**2+lam**2/2+2*s2*lam*L-(1+5*s2/2)*L**2-5*s2*polylog(2,s2-1)
print("Thm M1 - Prop(1/sqrt2):", nstr(abs(thm-M1_prop(1/s2)),3), " 8*M1 =", nstr(8*thm,22))
