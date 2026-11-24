from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, cos, sin, log, findroot
mp.dps = 60
j = 1j

def phi(s, N):   # exponent: (N-1) ln s - N ln(1+s) + i pi s/2
    return (N-1)*log(s) - N*log(1+s) + j*pi*s/2

def dphi(s, N):
    return (N-1)/s - N/(1+s) + j*pi/2

N = mpf(400)
theta = sqrt(pi*N)
rho = sqrt(2*N/pi)
omega = exp(j*pi/4)
# saddle from the quadratic: s^2 + (1+2i/pi) s - 2i(N-1)/pi = 0
disc = (1+2*j/pi)**2 + 8*j*(N-1)/pi
s_star = (-(1+2*j/pi) + disc**0.5)/2
print("s* =", s_star, " |s*| =", abs(s_star), "rho =", rho)
print("dphi(s*) =", dphi(s_star, N))
ps = phi(s_star, N)
print("phi(s*) =", ps)
print("-theta+i theta =", -theta, theta)
print("phi(s*) - (-theta+i theta) =", ps - (-theta + j*theta))
# s0 = rho*omega
s0 = rho*omega
print("phi(s0) =", phi(s0,N))
print("phi(s0) - (-theta+i theta) =", phi(s0,N) - (-theta+j*theta))
