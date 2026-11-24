from mpmath import mp, mpf, pi, sin, cos, quad, hyp2f1, gamma, sin as msin
import datetime

mp.dps = 50

def pi_p(p):
    p = mpf(p)
    return 2*pi/(p*sin(pi/p))

def F_p(y, p):
    """F(y) = int_0^y (1-t^p)^(-1/p) dt = y*2F1(1/p,1/p;1+1/p;y^p)"""
    p = mpf(p); y = mpf(y)
    return y*hyp2f1(1/p, 1/p, 1+1/p, y**p)

def Phi_p(y, p):
    """normalized: pi*F_p(y)/pi_p, maps [0,1]->[0,pi/2]"""
    return pi*F_p(y,p)/pi_p(p)

# Method 1: direct quadrature  b_n = (4/(n pi)) int_0^1 cos(n Phi_p(y)) dy
def b_n_quad(n, p):
    f = lambda y: cos(n*Phi_p(y,p))
    return 4/(n*pi)*quad(f, [0,1])

# Method 2: ODE integration. y(phi)=sin_p(phi): y'=(1-y^p)^(1/p), y(0)=0
# integrate to phi1=pi_p/2, alongside I' = cos_p(phi)*cos(n pi phi/pi_p)
from mpmath import odefun
def b_n_ode(n, p, tol=1e-40):
    p = mpf(p)
    phi1 = pi_p(p)/2
    def rhs(phi, Y):
        y, I = Y
        cosp = (1-y**p)**(1/p)
        return [cosp, cosp*cos(n*pi*phi/pi_p(p))]
    f = odefun(rhs, 0, [0,0], tol=tol, degree=40)
    Y1 = f(phi1)
    return 4/(n*pi)*Y1[1], Y1[0]

# sanity: p=2 -> delta_{n1}
for n in [1,3,5]:
    print("p=2 n=%d:"%n, b_n_quad(n,2))

# p=4, n=1,3: cross-check
for n in [1,3,5]:
    b1 = b_n_quad(n,4)
    b2, yend = b_n_ode(n,4)
    print(f"p=4 n={n}: quad={mp.nstr(b1,45)}")
    print(f"         ode ={mp.nstr(b2,45)}  yend-1={mp.nstr(yend-1,5)}")
