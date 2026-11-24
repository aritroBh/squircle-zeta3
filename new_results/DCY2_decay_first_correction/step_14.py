from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, cos, sin
mp.dps = 70
def Hc(s): return beta(s/2, mpf(3)/2)/2
omega = exp(1j*pi/4)
def I_tilde(N):
    f = lambda x: Hc(omega*x) * exp(1j*pi*omega*x/2) * (omega*x/(1+omega*x))**N / (omega*x) * omega
    rs = sqrt(2*N/pi)
    pts = [0, rs/8, rs/2, rs, 2*rs, 4*rs, 10*rs, mp.inf]
    return (2/pi) * mp.quad(f, pts)
varpi1 = ((24+12*pi+pi**2) + 1j*(24-12*pi+pi**2))/(48*sqrt(pi))

# cross-check at N=401 (k=400): certified 1.05786062767249754275690099276826282642433733510697906378343e-18
print("N=401 check:", mp.nstr(I_tilde(401).imag, 30))

for N in [100, 200, 400, 800, 1600, 3200]:
    It = I_tilde(N)
    M_main = It.imag
    theta = sqrt(pi*N)
    lead = -sqrt(pi)/N * exp(-theta)
    approx1 = lead * (cos(theta) + (varpi1.real*cos(theta)-varpi1.imag*sin(theta))/sqrt(N))
    res1 = (M_main - approx1) * N**2 * exp(theta)
    print(f"N={N:6d} res1*N^2*e^theta = {mp.nstr(res1,10)}")
