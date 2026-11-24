from mpmath import mp, mpf, mpc, pi, sqrt, exp, gamma, log, beta, cos, sin
mp.dps = 60

def Hc(s):   # CX6's H: (1/2) Beta(s/2, 3/2)
    return beta(s/2, mpf(3)/2)/2

omega = exp(1j*pi/4)

def I_tilde(N):
    # (2/pi) int_ray H(s) e^{i pi s /2} (s/(1+s))^N ds/s ; s = omega x
    f = lambda x: Hc(omega*x) * exp(1j*pi*omega*x/2) * (omega*x/(1+omega*x))**N / (omega*x) * omega
    rs = sqrt(2*N/pi)
    pts = [0, rs/8, rs/2, rs, 2*rs, 4*rs, 10*rs, mp.inf]
    return (2/pi) * mp.quad(f, pts)

# paper's varpi_1 closed form
varpi1 = ((24+12*pi+pi**2) + 1j*(24-12*pi+pi**2))/(48*sqrt(pi))
print("varpi1 =", mp.nstr(varpi1, 30))

for N in [100, 200, 400]:
    It = I_tilde(N)
    M_main = It.imag
    theta = sqrt(pi*N)
    lead = -sqrt(pi)/N * exp(-theta)
    approx0 = lead * cos(theta)
    approx1 = lead * (cos(theta) + (varpi1.real*cos(theta)-varpi1.imag*sin(theta))/sqrt(N))
    res0 = (M_main - approx0) * N**1.5 * exp(theta)
    res1 = (M_main - approx1) * N**2 * exp(theta)
    print(f"N={N} M_main={mp.nstr(M_main,12)}")
    print(f"   ratio M_main/lead0 = {mp.nstr(M_main/approx0,8)}   res0*N^1.5*e^theta = {mp.nstr(res0,8)}")
    print(f"   res1*N^2*e^theta = {mp.nstr(res1,8)}")
