from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j
def Hc(s): return beta(s/2, mpf(3)/2)/2

def Jnum_of(N):
    y = sqrt(pi/(2*N)); rho = 1/y
    omega = exp(j*pi/4); wb = exp(-j*pi/4)
    s0 = rho*omega
    def phi(s): return (N-1)*log(s) - N*log(1+s) + j*pi*s/2
    tmax = 1/sqrt(y)
    def integr(tau):
        t = tau/sqrt(y)
        return exp(phi(s0+omega*t)-phi(s0)) * Hc(s0+omega*t)/Hc(s0)
    return mp.quad(integr, [-tmax+mpf('0.0001'), 0, 5, mp.inf]), y

b = pi*exp(-j*pi/4)/2
for N in [mpf(400), mpf(2500), mpf(10000)]:
    Jn, y = Jnum_of(N)
    ratio = Jn/sqrt(pi/b)
    c1_est = (ratio-1)/y
    print(f"N={N}: J/sqrt(pi/b)-1 = {ratio-1},  c1_est = {c1_est}")
print("my c1 =", (1+1j)*(15-pi**2)/(8*sqrt(2)*pi))
