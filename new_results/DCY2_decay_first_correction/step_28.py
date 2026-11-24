from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j
def Hc(s): return beta(s/2, mpf(3)/2)/2

def model(N):
    y = sqrt(pi/(2*N)); rho = 1/y
    omega = exp(j*pi/4); wb = exp(-j*pi/4)
    s0 = rho*omega
    theta = sqrt(pi*N)
    phi_s0 = (N-1)*log(s0) - N*log(1+s0) + j*pi*s0/2
    # ell_m coefficients of ell(t) = (N-1)ln(1+yt) - N ln(1+yt/(1+y wb)) + i pi omega t/2
    def ell(m):
        if m == 1:
            return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
        return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
    l1, l2, l3, l4, l5 = ell(1), ell(2), ell(3), ell(4), ell(5)
    # h_j = (-omega)^j <u^j>_c / j!  with complex moments via Watson: exact via H derivatives
    # H^{(m)}(s) = (-1)^m int u^m sqrt(1-e^{-2u}) e^{-su} du -> compute by mpmath quadrature
    def Hder(m):
        f = lambda u: u**m * sqrt(1-exp(-2*u)) * exp(-s0*u)
        return (-1)**m * mp.quad(f, [0, mp.inf])
    H0 = Hc(s0)
    h = [1]
    for m in range(1,4):
        h.append(Hder(m)/H0 / mp.factorial(m))
    h1,h2,h3 = h[1],h[2],h[3]
    # moments mu_p^inf
    t0 = -l1/(2*l2)
    def mu(p):
        base = exp(-l2*t0**2) * sqrt(-pi/l2)
        tot = 0
        for m in range(0, p//2+1):
            tot += mp.binomial(p, 2*m) * t0**(p-2*m) * mp.factorial2(2*m-1) / (-2*l2)**m
        return base * tot
    # polynomial coefficients
    P = {}
    def add(p, c): P[p] = P.get(p, 0) + c
    add(0,1); add(1,h1); add(2,h2); add(3,h3)
    add(3,l3); add(4,l4); add(5,l5)
    add(6,l3**2/2); add(7,l3*l4); add(9,l3**3/6)
    add(4,h1*l3); add(5,h1*l4); add(5,h2*l3); add(7,h1*l3**2/2)
    M = sum(cp*mu(p) for p,cp in P.items())
    return (2/pi)*exp(phi_s0)*H0*omega*M

def I_tilde(N):
    f = lambda x: Hc(omega*x) * exp(j*pi*omega*x/2) * (omega*x/(1+omega*x))**N / (omega*x) * omega
    rs = sqrt(2*N/pi)
    return (2/pi) * mp.quad(f, [0, rs/8, rs/2, rs, 2*rs, 4*rs, 10*rs, mp.inf])

varpi1 = ((24+12*pi+pi**2) + 1j*(24-12*pi+pi**2))/(48*sqrt(pi))
for N in [100, 400, 2500]:
    N = mpf(N)
    theta = sqrt(pi*N)
    T1 = -1j*sqrt(pi)/N*exp(-theta+j*theta)*(1+varpi1/sqrt(N))
    It = I_tilde(N); Mo = model(N)
    print(f"N={int(N)}")
    print("  |It-model|  =", abs(It-Mo), " rel to N^-2 e^-th:", abs(It-Mo)*N**2*exp(theta))
    print("  |model - T1|=", abs(Mo-T1), " rel:", abs(Mo-T1)*N**2*exp(theta))
    print("  |It - T1|   =", abs(It-T1), " rel:", abs(It-T1)*N**2*exp(theta))
