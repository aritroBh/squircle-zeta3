import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

Nv = 10**10
N = mp.mpf(Nv); th = mp.sqrt(pi*N)
y = mp.sqrt(pi/(2*N))
f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
x0 = 1/(rho*y); t=rho*y
Phi0 = Phi(x0,N); p2 = Phi2(x0,N); g0 = g(x0)
om = mp.im(Phi1(x0,N))
c2 = abs(mp.re(p2))/2
c_w = mp.mpf('1.0')
w = mp.sqrt(mp.log(N)+c_w)
d = min(w/mp.sqrt(2*c2), x0/3)
sc = N**mp.mpf('-1.5')*mp.e**(-th)
P0 = abs(g0)*mp.e**mp.re(Phi0)
den = abs(p2)*d - om
first_true = 2*P0*mp.e**(-c2*d*d)/den
second_true = 2*P0*mp.e**(-c2*d*d)*abs(p2)/(den**2*2*c2*d)
print("TRUE pieces/sc: first =", mp.nstr(first_true/sc,12), " second =", mp.nstr(second_true/sc,12))

# code's certified first term: (2/pi)(pi/2)^1.5 * 2 rU^2.5 P0reg ecw / denreg_lo
G0reg = sqpi/4*2**mp.mpf('1.5')
eRe0_up = mp.e**((pi*rho**3/4)*(sqrt2/3)*y)   # mSreg_up*y approx
P0reg = G0reg*eRe0_up
ecw = mp.sqrt(2/pi)*mp.e**(-c_w/2)
p2reg = (pi/2)*rho**3*mp.sqrt(t**2+2*sqrt2*t+4)/(t**2+sqrt2*t+1)
c2reg = (pi*rho**3/4)*(t**3+2*sqrt2*t**2+4*t+sqrt2)/(1+sqrt2*t+t**2)**2
denreg = p2reg*w/mp.sqrt(2*c2reg) - (pi/(2*sqrt2))*(sqrt2*rho+rho**2*y-2*y*((rho-1)/y**2)-y**3*((rho-1)/y**2)**2)/(1+sqrt2*t+t**2)*mp.sqrt(y)
TMreg_first = 2*rho**mp.mpf('2.5')*P0reg*ecw/denreg
cert_first = (2/pi)*(pi/2)**mp.mpf('1.5')*TMreg_first
print("CERT pieces/sc: first =", mp.nstr(cert_first,12))
print("ratios: true_first/cert_first =", mp.nstr(first_true/sc/cert_first,12))
# dissect: factor-by-factor
print("P0*N^1.5 e^th =", mp.nstr(P0*N**mp.mpf('1.5')*mp.e**th,12), " vs t^2.5 P0reg (pi/2)^1.5 y^-0.5 =", mp.nstr(t**mp.mpf('2.5')*P0reg*(pi/2)**mp.mpf('1.5')/mp.sqrt(y),12))
print("e^{-c2 d^2} =", mp.nstr(mp.e**(-c2*d*d),12), " vs y*ecw =", mp.nstr(y*ecw,12))
print("den =", mp.nstr(den,12), " vs sqrt(y)*denreg =", mp.nstr(mp.sqrt(y)*denreg,12))
