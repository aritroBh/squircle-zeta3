from mpmath import mp, mpf, mpc, pi, sqrt, exp, gamma, log, im, re as mpre, cos, sin
mp.dps = 50

def Hf(s):
    return sqrt(pi)/2 * gamma(1+s/2)/gamma(mpf(3)/2+s/2)

omega = exp(1j*pi/4)

def I_tilde(N):
    # (2/pi) * int_ray H(s) E(s) s^{N-1}/(1+s)^N ds, s = omega*x, ds = omega dx
    f = lambda x: Hf(omega*x) * exp(1j*pi*omega*x/2) * (omega*x)**(N-1) / (1+omega*x)**N * omega
    return (2/pi) * mp.quad(f, [0, mp.inf])

for N in [100, 400]:
    It = I_tilde(N)
    theta = sqrt(pi*N)
    target = -1j*sqrt(pi)/N * exp(-theta+1j*theta)
    print("N=",N)
    print("  Itilde      =", It)
    print("  target lead =", target)
    print("  ratio       =", It/target)
