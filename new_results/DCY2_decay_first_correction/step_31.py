from mpmath import mp, mpf, pi, sqrt, exp, gamma, beta, log
mp.dps = 60
j = 1j
omega = exp(j*pi/4)
def df2(k):
    r = mpf(1)
    while k > 1: r *= k; k -= 2
    return r
def Hc(s): return beta(s/2, mpf(3)/2)/2

def pieces(N):
    y = sqrt(pi/(2*N)); rho = 1/y
    wb = exp(-j*pi/4)
    s0 = rho*omega
    def ell(m):
        if m == 1: return (N-1)*y - N*y/(1+y*wb) + j*pi*omega/2
        return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
    l1,l2,l3,l4,l5 = (ell(m) for m in range(1,6))
    def Hder(m):
        return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
    H0 = Hc(s0)
    h1,h2,h3 = (Hder(m)/H0/mp.factorial(m) for m in range(1,4))
    t0 = -l1/(2*l2)
    def mu(p):
        base = exp(-l2*t0**2) * sqrt(-pi/l2)
        return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
    terms = {
        '1': (0,1), 'h1': (1,h1), 'h2': (2,h2), 'h3': (3,h3),
        'l3': (3,l3), 'l4': (4,l4), 'l5': (5,l5),
        'l3^2/2': (6,l3**2/2), 'l3l4': (7,l3*l4), 'l3^3/6': (9,l3**3/6),
        'h1l3': (4,h1*l3), 'h1l4': (5,h1*l4), 'h2l3': (5,h2*l3), 'h1l3^2/2': (7,h1*l3**2/2),
    }
    out = {}
    for name,(p,c) in terms.items():
        out[name] = c*mu(p)
    return out, y

N = mpf(400)
out, y = pieces(N)
mu0 = None
tot = sum(out.values())
print(f"N={int(N)}, y={y}")
for name, v in out.items():
    print(f"  {name:10s} {v}   |v|/y = {abs(v)/y}")
