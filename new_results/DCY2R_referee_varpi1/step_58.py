import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp
mp.dps = 90
j=1j; omega = exp(j*pi/4); SQ2 = sqrt(2)

def Hc(s):
    # stable: sqrt(pi)/4 Gamma(s/2)/Gamma(3/2+s/2)
    return sqrt(pi)/4*mp.gamma(s/2)/mp.gamma(mpf(3)/2+s/2)

N=100000
N=mpf(N); y = sqrt(pi/(2*N)); rho=1/y; s0 = rho*omega
kap2 = mpf('2.48')**2 + 4*mp.log(N/400); kap = sqrt(kap2)
c2m = mpf('0.9')*(pi/(2*SQ2))*y; D = kap/sqrt(2*c2m)
H0 = Hc(s0)
def Hder(m):
    return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
h = [omega**m*Hder(m)/H0/mp.factorial(m) for m in range(1,4)]
H_lo = sqrt(pi/2)*y**1.5*(1-3*y/4-y**2)
x = (rho-D)/SQ2
M4 = SQ2*mp.gamma(mpf(11)/2)*x**(-mpf(11)/2)/H_lo
worst = 0
for k in range(41):
    t = -D + 2*D*k/40
    F = Hc(s0+omega*t)/H0
    PH = 1 + h[0]*t + h[1]*t**2 + h[2]*t**3
    r = abs(F-PH)/(M4*abs(t)**4/24) if t != 0 else 0
    worst = max(worst, r)
print(f"N=100000: max ratio = {mp.nstr(worst,6)} (dps=90)")
