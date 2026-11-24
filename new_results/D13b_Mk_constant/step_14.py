
import mpmath as mp
mp.mp.dps = 60

def B(z):
    return 1 - mp.exp(z)/(1+mp.sqrt(1+mp.exp(2*z)))

def Disc(s, r, dlt=mp.mpf('1e-40')):
    zs = s*1j*mp.pi/2
    return B(zs + r + 1j*dlt) - B(zs + r - 1j*dlt)

N = mp.mpf(401)
s = 1
zs = s*1j*mp.pi/2

# exact saddle of phi(u) = -u - N log(1+zs/u):  u(u+zs) = N zs
us = (-zs + mp.sqrt(zs**2 + 4*N*zs))/2
print("u_s =", mp.nstr(us, 20))

# ---------- STEP A: zeta-integral at u = u_s ----------
def Jhat_exact(u, R=mp.mpf('2.5')):
    f = lambda r: Disc(s, r) * (1+(zs+r)/u)**(-N)
    return mp.quad(f, [0, R])/(2*mp.pi*1j)

lam = N/(us+zs)
c3 = mp.mpf('-0.5'); c5 = mp.mpf(5)/24
Jhat_asymp = (mp.e**(-N*mp.log(1+zs/us)) * mp.gamma(mp.mpf('1.5')) * lam**(-mp.mpf('1.5'))
              * (1 + mp.mpf('1.5')*c3/lam + mp.mpf('15')/8/N + c5*mp.mpf('15')/4/lam**2))
Je = Jhat_exact(us)
print("J_exact(u_s) =", mp.nstr(Je, 18))
print("J_asymp(u_s) =", mp.nstr(Jhat_asymp, 18))
print("ratio exact/asymp =", mp.nstr(Je/Jhat_asymp, 14))
print("expected ratio = 1 + O(N^-3/2) ~ 1 +-", mp.nstr(mp.mpf(1)/N**mp.mpf('1.5'),3))
