import mpmath as mp, time
mp.mp.dps = 60
beta = mp.pi**2/8; m = mp.mpf(2); cm = 1-1/m; ym = cm*(1-mp.log(cm))
def Uroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, -1)
def Vroot(y): return -cm*mp.lambertw(-mp.exp(-y/cm)/cm, 0)
def ell(y):
    if y <= ym: return mp.mpf(0)
    return Uroot(y)**(1/m) - Vroot(y)**(1/m)
def phi_s(s,c):
    r = mp.sqrt(1+2*c*s+s**2)
    return (1-c**2)/(r*(r+s+c))
th = mp.pi/3; c = mp.cos(th); n = 10**6
Nn = mp.sqrt(beta*mp.mpf(n)); Lam = mp.log(m*Nn)
ustar = Lam - ym
print("Lam =", mp.nstr(Lam,12), " singularity at u* =", mp.nstr(ustar,12))
# Route A (mine): split at 0
i1 = mp.quad(lambda u: phi_s(mp.exp(-u),c)*mp.exp(-u)*ell(Lam+u), [0, mp.inf])
i2 = mp.quad(lambda u: phi_s(mp.exp(u),c)*mp.exp(u)*ell(Lam-u), [0, mp.inf])
IA = i1+i2
# Route B (theirs): [-60..60] with their breakpoints
IB = mp.quad(lambda u: ell(Lam-u)*phi_s(mp.exp(u),c)*mp.exp(u), [-60, Lam-12, Lam-4, Lam, Lam+4, 60], maxdegree=10)
# Route C: careful, breakpoint AT the singularity, u = log s
IC = mp.quad(lambda u: ell(Lam-u)*phi_s(mp.exp(u),c)*mp.exp(u), [-80, 0, ustar-4, ustar, ustar+2, ustar+8, 60], maxdegree=12)
# Route D: direct s-integral with split at s0 = exp(-(Lam-ym))
s0 = mp.exp(-(Lam-ym))
ID = mp.quad(lambda s: phi_s(s,c)*ell(Lam+mp.log(1/s)), [0, s0/10, s0, 10*s0, mp.inf], maxdegree=12)
print("A:", mp.nstr(IA,25))
print("B:", mp.nstr(IB,25))
print("C:", mp.nstr(IC,25))
print("D:", mp.nstr(ID,25))
print("B-A:", mp.nstr(IB-IA,3), " C-A:", mp.nstr(IC-IA,3), " D-A:", mp.nstr(ID-IA,3))
