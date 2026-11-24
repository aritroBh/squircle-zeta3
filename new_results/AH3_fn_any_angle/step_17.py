
import mpmath as mp, time
mp.mp.dps = 55

def Lambda_int(nu, th):  # noninteger order
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)
def Lambda_poly(m, th):  # integer order polynomial in c
    c = mp.cos(th)
    return (mp.legendre(m-1, c) - mp.legendre(m+1, c))/(2*m+1)

def exact_deficit_x(th, n, M=45):
    # Thm P2: L = 1+c + (pi eps/sin pi eps) Lam_eps n^{-eps} + sum_m (-1)^m eps/(m-eps) Lam_m n^{-m}
    nn = mp.mpf(n); c = mp.cos(th); eps = 1/(nn-1)
    L = 1 + c + (mp.pi*eps/mp.sin(mp.pi*eps))*Lambda_int(eps, th)*nn**(-eps)
    for m in range(1, M+1):
        L += (-1)**m * eps/(m-eps) * Lambda_poly(m, th) * nn**(-m)
    return 2 - L

def hc(s, th):
    c = mp.cos(th)
    return 1 + s - mp.sqrt(1 + 2*c*s + s**2)

def quad_deficit_x(th, n):
    nn = mp.mpf(n)
    def integrand(v):
        t = mp.exp(v)/nn
        yp = nn*mp.expm1((nn-1)*mp.log1p(-t)) + 0  # careful: use exp directly
        yp = nn*mp.exp((nn-1)*mp.log1p(-t))
        return hc(yp, th)*mp.exp(v)/nn
    return mp.quad(integrand, [-130, 0, 3, mp.log(nn)], maxdegree=8)

t0=time.time()
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3]:
    for n in [100, 1000, 10000]:
        q = quad_deficit_x(th, n); e = exact_deficit_x(th, n)
        print(f" th={mp.nstr(th,6)} n={n}: quad={mp.nstr(q,18)} exact={mp.nstr(e,18)} |diff|={mp.nstr(abs(q-e),3)}")
print('elapsed', time.time()-t0)
t0=time.time()
q = quad_deficit_x(mp.pi/2, 100000); e = exact_deficit_x(mp.pi/2, 100000)
print('n=1e5:', mp.nstr(q,18), mp.nstr(e,18), 'diff', mp.nstr(abs(q-e),3), 'time', time.time()-t0)
