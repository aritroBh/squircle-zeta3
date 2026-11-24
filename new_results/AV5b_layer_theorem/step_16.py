
import mpmath as mp
mp.mp.dps = 80

def S(phi, m, p):
    return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
def rho(phi, m, p):
    return S(phi, m, p)**(-mp.mpf(1)/p)
def speed(phi, m, p):
    r = rho(phi, m, p)
    dr = mp.diff(lambda ph: rho(ph, m, p), phi)
    return mp.sqrt(r*r + dr*dr)

# m=3, p=1: one edge sector: between adjacent vertex directions pi/2 and pi/2-pi/3 = pi/6
# edge midpoint at pi/3. Exact edge length: 2*tan(pi/6)*sin(pi/6)=... rho_v=tan(pi/6)=1/sqrt3, L=2*rho_v*sin(pi/6)=1/sqrt3
m = 3; p = mp.mpf(1)
seg, err = mp.quad(lambda ph: speed(ph, m, p), [mp.pi/6, mp.pi/3, mp.pi/2], error=True)
print("one edge sector integral:", mp.nstr(seg, 30), "err est", mp.nstr(err, 5))
print("exact edge length 1/sqrt3:", mp.nstr(1/mp.sqrt(3), 30))

# try finer quadrature control
mp.mp.dps = 60
seg2 = mp.quad(lambda ph: speed(ph, m, p), [mp.pi/6, mp.pi/3, mp.pi/2], method="gauss-legendre")
print("gauss-legendre:", mp.nstr(seg2, 30))
