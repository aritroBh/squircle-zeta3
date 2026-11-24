from mpmath import mp, mpf, pi, sin, cos, quad, sqrt, odefun
import time
mp.dps = 40

def pi_p(p):
    p = mpf(p)
    return 2*pi/(p*sin(pi/p))

def F_p(y, p):
    p=mpf(p); y=mpf(y)
    if y == 0: return mpf(0)
    if y == 1: return pi_p(p)/2
    return y*quad(lambda s: (1-(y*s)**p)**(-1/p), [0,1])

def Phi_p(y,p):
    return pi*F_p(y,p)/pi_p(p)

t0=time.time()
print("Phi test:", mp.nstr(Phi_p(0.5,4),30), f"{time.time()-t0:.2f}s")
t0=time.time()
print("Phi near 1:", mp.nstr(Phi_p(0.999999,4),30), f"{time.time()-t0:.2f}s")

# Method 1
def b_n_quad(n, p):
    return 4/(n*pi)*quad(lambda y: cos(n*Phi_p(y,p)), [0,1])

# Method 2: psi-parametrization, G(psi)=(2/p) int_0^psi sin^a u cos^b u du
def b_n_psi(n,p):
    p=mpf(p); a=2/p-1; b=1-2/p
    pp = pi_p(p)
    def G(psi):
        if psi==0: return mpf(0)
        return (2/p)*quad(lambda u: sin(u)**a * cos(u)**b, [0,psi])
    return 8/(n*p*pi)*quad(lambda psi: sin(psi)**a*cos(psi)*cos(n*pi*G(psi)/pp), [0,pi/2])

t0=time.time(); b1 = b_n_quad(1,4); print("m1 b1(p=4):", mp.nstr(b1,35), f"{time.time()-t0:.1f}s")
t0=time.time(); b2 = b_n_psi(1,4);  print("m2 b1(p=4):", mp.nstr(b2,35), f"{time.time()-t0:.1f}s")
