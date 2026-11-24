# step_20.py : z = x^n y^{ρn}: exact identity E = n J_ρ - C_ρ ; numerical validation
from mpmath import mp, mpf, sqrt, log, nstr, pi, zeta, gamma, inf, quad, exp, cosh, binomial
import math
mp.dps = 60
L2 = log(2); z2 = pi**2/6; z3 = zeta(3); z4 = pi**4/90

def J_series(n, rho, kmax=200):
    n = mpf(n); m = rho*n
    s = mpf(0)
    for k in range(kmax):
        b = binomial(mpf('0.5'), k)
        t1 = rho**(2*k)/(m-2*k+1)*(1/(n+2*k) - rho**(m-2*k+1)/(m+n+1))
        t2 = rho**(m+1)/((m+2*k)*(m+n+1))
        s += b*(t1+t2)
    return s

def E_direct(n, rho):
    # A(n) = ∬ sqrt(1 + n² x^{2n-2} y^{2m} + m² x^{2n} y^{2m-2}) dxdy,  m=ρn ; E = A-1
    n = mpf(n); m = rho*n
    def F(x,y):
        return sqrt(1 + (n*x**(n-1)*y**m)**2 + (m*x**n*y**(m-1))**2) - 1
    # integrate via tanh-sinh on transformed coords x=u^{1/n}, y=v^{1/n}: smooth
    # use the planner substitution directly: E = (1/(ρn²))∬[√(1+W²)-1] e^{-eu}e^{-ew/ρ} dudw
    e = 1/n
    def G(u,w):
        W = n*exp(-(u+w))*sqrt(exp(2*e*u) + rho**2*exp(2*e*w/rho))
        return (sqrt(1+W*W)-1)*exp(-e*u)*exp(-e*w/rho)
    L = log(n)
    I = quad(G, [0, L/2, L, L+30], [0, L/2, L, L+30])
    return I/(rho*n**2)

def fW(w):
    w = mpf(w); return 1 + w - sqrt(1+w*w)
def C_rho_int(n, rho):
    # C_ρ = (1/(ρn²))∬ f(W) e^{-eu}e^{-ew/ρ} dudw
    n = mpf(n); e = 1/n
    def G(u,w):
        W = n*exp(-(u+w))*sqrt(exp(2*e*u) + rho**2*exp(2*e*w/rho))
        return fW(W)*exp(-e*u)*exp(-e*w/rho)
    L = log(n)
    I = quad(G, [0, L/2, L, L+30], [0, L/2, L, L+30])
    return I/(rho*n**2)

for n, rho in [(100, mpf('0.5')), (100, mpf('2')), (1000, mpf('0.5'))]:
    Ed = E_direct(n, rho)
    Eid = n*J_series(n, rho) - C_rho_int(n, rho)
    print("n=%s rho=%s: E direct = %s" % (n, rho, nstr(Ed, 30)))
    print("          E = nJ - C_ρ = %s   diff = %s" % (nstr(Eid, 30), nstr(Ed-Eid, 8)))
