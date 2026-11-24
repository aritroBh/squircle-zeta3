
from mpmath import mp, mpf, pi, sin, cos, tan, log, quad, sqrt, sin as msin

def pi_p(p):
    p = mpf(p); return 2*pi/(p*sin(pi/p))

# q(p) = pi/pi_p = (p/2) sin(pi/p); derivatives at p=2:
Q0 = mpf(1)
Q1 = mpf('0.5')
Q2 = -pi**2/16
Q3 = 3*pi**2/32

def _gs(y, c2):
    s2 = 1-c2; d = -c2
    if y < 1:
        ls = log(y)
    else:
        ls = -c2/2 - c2*c2/8 - c2**3/24
    lc2 = log(c2)
    g1 = (-2*s2*ls + d*lc2)/(4*d)
    g2 = (2*s2*s2*ls**2 - 2*s2*d*(ls-1)*ls - d*d*lc2)/(4*d*d)
    g3 = (-8*s2**3*ls**3 + 6*s2*s2*d*(2*ls-1)*ls**2 + 6*s2*d*d*(ls-1)*ls
          - 4*s2*d*d*ls**3 + 3*d**3*lc2)/(8*d**3)
    return g1, g2, g3

def H1u(u):
    s = sin(u); c = cos(u); c2 = c*c
    if s == 0: return mpf(0)
    g1,g2,g3 = _gs(s,c2)
    H = c2**mpf('-0.5')
    return H*g1*c

def H2u(u):
    s = sin(u); c = cos(u); c2 = c*c
    if s == 0: return mpf(0)
    g1,g2,g3 = _gs(s,c2)
    H = c2**mpf('-0.5')
    return H*(g2+g1*g1)*c

def H3u(u):
    s = sin(u); c = cos(u); c2 = c*c
    if s == 0: return mpf(0)
    g1,g2,g3 = _gs(s,c2)
    H = c2**mpf('-0.5')
    return H*(g3+3*g1*g2+g1**3)*c

def L1(phi): return quad(H1u, [0, phi])
def L2(phi): return quad(H2u, [0, phi])
def L3(phi): return quad(H3u, [0, phi])

def Psi1(phi): return Q0*L1(phi) + Q1*phi
def Psi2(phi): return Q0*L2(phi) + 2*Q1*L1(phi) + Q2*phi
def Psi3(phi): return Q0*L3(phi) + 3*Q1*L2(phi) + 3*Q2*L1(phi) + Q3*phi

def make_psis(need3=True, need2=True):
    cache = {}
    def Ls(phi):
        if phi not in cache:
            l1 = quad(H1u,[0,phi])
            l2 = quad(H2u,[0,phi]) if (need2 or need3) else None
            l3 = quad(H3u,[0,phi]) if need3 else None
            cache[phi] = (l1,l2,l3)
        return cache[phi]
    def P1(phi):
        l1,l2,l3 = Ls(phi); return Q0*l1 + Q1*phi
    def P2(phi):
        l1,l2,l3 = Ls(phi); return Q0*l2 + 2*Q1*l1 + Q2*phi
    def P3(phi):
        l1,l2,l3 = Ls(phi); return Q0*l3 + 3*Q1*l2 + 3*Q2*l1 + Q3*phi
    return P1, P2, P3

def bn_der(n, k):
    n = mpf(n)
    P1, P2, P3 = make_psis(need3=(k==3), need2=(k>=2))
    if k == 1:
        f = lambda u: -n*sin(n*u)*P1(u)*cos(u)
    elif k == 2:
        f = lambda u: (-n**2*cos(n*u)*P1(u)**2 - n*sin(n*u)*P2(u))*cos(u)
    elif k == 3:
        f = lambda u: (n**3*sin(n*u)*P1(u)**3 - 3*n**2*cos(n*u)*P1(u)*P2(u)
                       - n*sin(n*u)*P3(u))*cos(u)
    # split: [0, a] with panels, [a, pi/2] with u = pi/2 - t^2 transform
    a = mpf('1.2')
    pts = [j*a/int(n) for j in range(int(n)+1)]
    I1 = quad(f, pts)
    b = sqrt(pi/2 - a)
    g = lambda t: f(pi/2 - t*t)*2*t
    I2 = quad(g, [0, b])
    return 4/(n*pi)*(I1+I2)

def F_p(y, p):
    p=mpf(p); y=mpf(y)
    if y == 0: return mpf(0)
    if y == 1: return pi_p(p)/2
    return y*quad(lambda s: (1-(y*s)**p)**(-1/p), [0,1])

def Phi_p(y,p):
    return pi*F_p(y,p)/pi_p(p)

def b_n_quad(n, p):
    return 4/(n*pi)*quad(lambda y: cos(n*Phi_p(y,p)), [0,1])
