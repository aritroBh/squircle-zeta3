
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
mp.dps = 120
u4 = exp(mpc(0,1)*mp.pi/4)
s2 = mp.sqrt(2)

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def Jm(N, Z, m, p):
    a = abs(Z); sq2 = mp.sqrt(2)
    f0 = lambda t: p/sq2 - (N/2)*(sq2*a*t + 2*a*a)/(t**3 + sq2*a*t*t + a*a*t)
    ts = mp.findroot(f0, mp.sqrt(N*a/p))
    Rp = ts*p/sq2 + (N/2)*mp.log(1 + sq2*a/ts + (a/ts)**2)
    esc = exp(Rp)
    def ft(t):
        x = t*u4
        return exp(-p*x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x * u4 * esc
    w = ts/6
    pts = [0, ts*mpf('0.3'), ts*mpf('0.6'), ts-w, ts-w/3, ts, ts+w/3, ts+w, ts*mpf('1.6'), ts*mpf('3'), ts*mpf('10'), mp.inf]
    return mp.quadgl(ft, pts)/esc

def Iplus_m(N, a, p, m):
    Z = 1j*a
    bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    return (+1j/mp.pi)*gam(a)*gms[m]*bf*Jm(N,Z,m,p)

def Iplus(N, a, p):
    return sum(Iplus_m(N,a,p,m) for m in range(3))

# model pieces
def y_p(sig, p, q):
    sp_ = mp.sqrt(p)*sig
    return mp.sqrt(1+q*q*sp_**2/4) - q*sp_/2
def x_p(sig, p, q):
    sp_ = mp.sqrt(p)*sig
    return q*y_p(sig,p,q)/sp_
def phi_p(sig, p, q, Z):
    xp = x_p(sig,p,q)
    return p*xp + sig**-2*mp.log(1+Z/xp)
def phi2der(sig, p, q, Z):
    xp = x_p(sig,p,q)
    return sig**-2*(1/xp**2 - 1/(xp+Z)**2)
def B0(s, q): return 1 + s/(2*q) - mp.mpf(5)/8*s**2
def B1(s, q): return 1 + 3*s/(2*q) + (-mp.mpf(7)/4 + 3/(4*q**2))*s**2

a = mp.pi/2; Z = 1j*a; q = mp.sqrt(Z); g1 = g1f(a)
for N in [800, 2000, 5000, 12000]:
    sig = 1/mp.sqrt(mpf(N))
    th = mp.sqrt(2*a*N)
    I1 = Iplus(N, a, 1); I2 = Iplus(N, a, 2)
    R_quad = -2*I2/I1
    # model
    lead = 2*(s2-1)*q/sig
    D = phi_p(sig,2,q,Z) - phi_p(sig,1,q,Z)
    S_exp = exp(-D)
    w1 = x_p(sig,1,q)/y_p(sig,1,q)**2
    w2 = x_p(sig,2,q)/y_p(sig,2,q)**2
    Pref = sqrt(phi2der(sig,1,q,Z)/phi2der(sig,2,q,Z))*(x_p(sig,1,q)/x_p(sig,2,q))*(w2/w1)**mp.mpf('1.5')
    F3 = B0(s2*sig,q)/B0(sig,q)
    Nf = sig**-2
    m1_1 = mp.mpf('1.5')*g1*w1/Nf * B1(sig,q)/B0(sig,q)
    m1_2 = mp.mpf('1.5')*g1*w2/Nf * B1(s2*sig,q)/B0(s2*sig,q)
    F4 = (1+m1_2)/(1+m1_1)
    R_model = -2*S_exp*Pref*F3*F4
    print('N=%5d sig=%.4f' % (N, float(sig)))
    print('  R_quad *e^lead =', mp.nstr(R_quad*exp(lead), 16))
    print('  R_model*e^lead =', mp.nstr(R_model*exp(lead), 16))
