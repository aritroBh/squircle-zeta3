
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma
mp.dps = 120
u4 = exp(mpc(0,1)*mp.pi/4)

def gam(a): return mp.sqrt(2*mp.sin(a))*mp.e**(1j*(mp.pi/4 - a/2))
def g1f(a): return -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
def g2f(a):
    s = mp.sin(a); E = mp.e**(1j*a)
    return mp.mpf(29)/96 + E**2*(mp.mpf(1)/2 + 1/(32*s**2)) + (1j*E/s)*(E**2/4 - mp.mpf(3)/16)

def Jm(N, Z, m, p):   # weight e^{-p x}
    a = abs(Z); s2 = mp.sqrt(2)
    # Re phi(t) = p t/sqrt2 + (N/2) log(1+sqrt2 a/t + a^2/t^2)
    f0 = lambda t: p/s2 - (N/2)*(s2*a*t + 2*a*a)/(t**3 + s2*a*t*t + a*a*t)
    ts = mp.findroot(f0, mp.sqrt(N*a/p))
    Rp = ts*p/s2 + (N/2)*mp.log(1 + s2*a/ts + (a/ts)**2)
    esc = exp(Rp)
    def ft(t):
        x = t*u4
        return exp(-p*x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x * u4 * esc
    w = ts/6
    pts = [0, ts*mpf('0.3'), ts*mpf('0.6'), ts-w, ts-w/3, ts, ts+w/3, ts+w, ts*mpf('1.6'), ts*mpf('3'), ts*mpf('10'), mp.inf]
    return mp.quadgl(ft, pts)/esc

def Iplus(N, a, p):
    Z = 1j*a; tot = 0
    gms = [mp.mpf(1), g1f(a), g2f(a)]
    for m in range(3):
        bf = mp.gamma(m+mp.mpf('1.5'))*mp.gamma(N-m-mp.mpf('1.5'))/mp.gamma(N)
        tot += gms[m]*bf*Jm(N,Z,m,p)
    return (+1j/mp.pi)*gam(a)*tot

for a in [mp.pi/2, 5*mp.pi/6]:
    print('alpha =', mp.nstr(a, 6))
    for N in [800, 2000, 5000]:
        th = mp.sqrt(2*a*N)
        I1 = Iplus(N, a, 1)
        I2 = Iplus(N, a, 2)          # this is the raw p=2 integral (without the weight coefficient -2)
        # image contribution = (-2)*I2 ; relative to main I1:
        R = -2*I2/I1
        C = R*exp((mp.sqrt(2)-1)*th*(1+1j))   # should -> constant ~ -(1+c sig+...)
        sig = 1/mp.sqrt(mpf(N))
        print('  N=%5d sig=%.4f  C_img = %s' % (N, float(sig), mp.nstr(C, 18)))
