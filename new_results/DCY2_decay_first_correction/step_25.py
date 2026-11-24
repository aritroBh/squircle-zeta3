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
# Richardson extrapolate c1(N)=(J/sqrt(pi/b)-1)/y ~ c1 + d1*y
vals = []
for N in [mpf(1e5), mpf(4e5), mpf(1.6e6)]:
    Jn, y = Jnum_of(N)
    vals.append((y, (Jn/sqrt(pi/b)-1)/y))
for yv, cv in vals: print(yv, cv)
# quadratic extrapolation in y
import numpy as np
A = np.array([[1,float(yv),float(yv)**2] for yv,_ in vals])
cs = np.array([complex(cv) for _,cv in vals])
sol = np.linalg.solve(A, cs)
print("extrapolated c1_true:", sol[0])

# pieces (symbolic) numerically
wb = exp(-j*pi/4)
def avg(f):  # <f(tau)> Gaussian
    num = mp.quad(lambda t: f(t)*exp(-b*t**2), [-mp.inf,0,mp.inf])
    den = sqrt(pi/b)
    return num/den
G1 = lambda t: -(1 + pi*wb**2/2)*t + (pi*wb/2)*t**3
G2 = lambda t: (mpf(1)/2 + 3*pi*wb**2/4)*t**2 - (pi*wb/2)*t**4
h1, h2 = mpf(-3)/2, mpf(15)/8
print("avg G2        =", avg(G2))
print("avg G1^2/2    =", avg(lambda t: G1(t)**2/2))
print("avg h2 t^2    =", avg(lambda t: h2*t**2))
print("avg G1 h1 t   =", avg(lambda t: G1(t)*h1*t))
tot = avg(G2)+avg(lambda t: G1(t)**2/2)+avg(lambda t: h2*t**2)+avg(lambda t: G1(t)*h1*t)
print("sum =", tot)
