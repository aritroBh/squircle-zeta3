
import mpmath as mp
mp.mp.dps = 60

def B(z):
    return 1 - mp.exp(z)/(1+mp.sqrt(1+mp.exp(2*z)))

def Disc(s, r, dlt=mp.mpf('1e-40')):
    zs = s*1j*mp.pi/2
    return B(zs + r + 1j*dlt) - B(zs + r - 1j*dlt)

N = mp.mpf(401); s = 1
zs = s*1j*mp.pi/2
us = (-zs + mp.sqrt(zs**2 + 4*N*zs))/2

# LHS: (1/2pi i) ∫ Disc(r) exp(-N log(1 + r/(us+zs))) dr   [exact]
f = lambda r: Disc(s, r)*mp.e**(-N*mp.log(1 + r/(us+zs)))
LHS = mp.quad(f, [0, mp.mpf('2.5')])/(2*mp.pi*1j)

# RHS: (1/2pi i) ∫ [2 s sqrt2 (zeta^{1/2} - 1/2 zeta^{3/2} + 5/24 zeta^{5/2} ...)] e^{-lam zeta} dzeta
lam = N/(us+zs)
RHS0 = mp.gamma(mp.mpf('1.5'))*lam**mp.mpf('-1.5')
RHS1 = -mp.mpf('0.5')*mp.gamma(mp.mpf('2.5'))*lam**mp.mpf('-2.5')
RHS2 = mp.mpf(5)/24*mp.gamma(mp.mpf('3.5'))*lam**mp.mpf('-3.5')
DiscL = 2*s*mp.sqrt(2)
RHS = DiscL*(RHS0+RHS1+RHS2)/(2*mp.pi*1j)*2*mp.pi*1j   # note: (1/2pii)∫ = we keep explicit below
RHS = DiscL*(RHS0+RHS1+RHS2)/(2*mp.pi*1j)
print("lam =", mp.nstr(lam, 14))
print("LHS (exact) =", mp.nstr(LHS, 16))
print("RHS (asymp) =", mp.nstr(RHS, 16))
print("ratio =", mp.nstr(LHS/RHS, 14))
print()
# also check Disc series against numeric at one r
r0 = mp.mpf('0.05')
ser = 2*s*mp.sqrt(2)*(mp.sqrt(r0) - mp.mpf('0.5')*r0**mp.mpf('1.5') + mp.mpf(5)/24*r0**mp.mpf('2.5'))
print("Disc numeric at r=0.05:", mp.nstr(Disc(s, r0), 16))
print("Disc series           :", mp.nstr(ser, 16))
