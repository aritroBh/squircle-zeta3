from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, cos, sin, log
mp.dps = 60
j = 1j
N = mpf(400)
y = sqrt(pi/(2*N))
rho = 1/y
omega = exp(j*pi/4); wb = exp(-j*pi/4)
s0 = rho*omega
theta = sqrt(pi*N)

def Hc(s): return beta(s/2, mpf(3)/2)/2
def phi(s):
    return (N-1)*log(s) - N*log(1+s) + j*pi*s/2

# exact tau-integral: s = s0 + omega*t, t = tau/sqrt(y)
def G_exact(tau):
    t = tau/sqrt(y)
    return phi(s0+omega*t) - phi(s0)
def Hratio(tau):
    t = tau/sqrt(y)
    return Hc(s0+omega*t)/Hc(s0)

b = pi*wb/2
Jnum = mp.quad(lambda tau: exp(G_exact(tau))*Hratio(tau), [-mp.inf, 0, mp.inf])
print("J exact =", Jnum)

# formal: sqrt(pi/b)*(1 + c1*y)
c1 = (1+1j)*(15-pi**2)/(8*sqrt(2)*pi)
Jformal = sqrt(pi/b)*(1 + c1*y)
print("J formal =", Jformal)
print("ratio =", Jnum/Jformal)

# prefactor check
P_exact = (2/pi)*exp(phi(s0))*Hc(s0)*omega/sqrt(y)
P0 = -1j*(2/sqrt(pi))*y**2*exp(-theta+j*theta)
S1 = pi*(1+1j)/(6*sqrt(2))
pref_corr = S1 - 3*wb/4
print("prefactor ratio =", P_exact/(P0*(1+pref_corr*y)))
print("P_exact/P0 =", P_exact/P0)
print("1+pref_corr*y =", 1+pref_corr*y)
