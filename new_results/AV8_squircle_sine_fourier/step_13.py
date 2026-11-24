module = r'''
from mpmath import mp, mpf, pi, sin, cos, tan, log, quad, sqrt, sin as msin

def pi_p(p):
    p = mpf(p); return 2*pi/(p*sin(pi/p))

# q(p) = pi/pi_p = (p/2) sin(pi/p); derivatives at p=2:
Q0 = mpf(1)
Q1 = mpf('0.5')
Q2 = -pi**2/16
Q3 = 3*pi**2/32

def H1u(u):
    s = sin(u); c = cos(u); s2 = s*s; c2 = 1-s2
    A = 2*s2*log(s) - (s2-1)*log(c2)
    return A/(4*c2**mpf('1.5'))*c

def H2u(u):
    s = sin(u); c = cos(u); s2 = s*s; c2 = 1-s2
    ls = log(s); lc2 = log(c2)
    A = 2*s2*ls - (s2-1)*lc2
    num = 8*s2*s2*ls**2 - 8*s2*(s2-1)*(ls-1)*ls - 4*c2**2*lc2 + A**2
    return num/(16*c2**mpf('2.5'))*c

def H3u(u):
    s = sin(u); c = cos(u); s2 = s*s; c2 = 1-s2
    ls = log(s); lc2 = log(c2)
    A = 2*s2*ls - (s2-1)*lc2
    B = -2*s2*s2*ls**2 + 2*s2*(s2-1)*(ls-1)*ls + c2**2*lc2
    num = (64*s2**3*ls**3 - 48*s2*s2*(s2-1)*(2*ls-1)*ls**2 - 48*s2*(s2-1)**2*(ls-1)*ls
           + 32*s2*(s2-1)**2*ls**3 - 24*c2**3*lc2 + A**3 - 12*A*B)
    return num/(64*c2**mpf('3.5'))*c

def L1(phi): return quad(H1u, [0, phi])
def L2(phi): return quad(H2u, [0, phi])
def L3(phi): return quad(H3u, [0, phi])

def Psi1(phi): return Q0*L1(phi) + Q1*phi
def Psi2(phi): return Q0*L2(phi) + 2*Q1*L1(phi) + Q2*phi
def Psi3(phi): return Q0*L3(phi) + 3*Q1*L2(phi) + 3*Q2*L1(phi) + Q3*phi

def bn_der(n, k):
    """b_n^{(k)}(2) for k=1,2,3, odd n"""
    n = mpf(n)
    if k == 1:
        f = lambda u: -n*sin(n*u)*Psi1(u)*cos(u)
    elif k == 2:
        f = lambda u: (-n**2*cos(n*u)*Psi1(u)**2 - n*sin(n*u)*Psi2(u))*cos(u)
    elif k == 3:
        f = lambda u: (n**3*sin(n*u)*Psi1(u)**3 - 3*n**2*cos(n*u)*Psi1(u)*Psi2(u)
                       - n*sin(n*u)*Psi3(u))*cos(u)
    return 4/(n*pi)*quad(f, [0, pi/2])

def F_p(y, p):
    p=mpf(p); y=mpf(y)
    if y == 0: return mpf(0)
    if y == 1: return pi_p(p)/2
    return y*quad(lambda s: (1-(y*s)**p)**(-1/p), [0,1])

def Phi_p(y,p):
    return pi*F_p(y,p)/pi_p(p)

def b_n_quad(n, p):
    return 4/(n*pi)*quad(lambda y: cos(n*Phi_p(y,p)), [0,1])
'''
open('squircle.py','w').write(module)
print("written")
