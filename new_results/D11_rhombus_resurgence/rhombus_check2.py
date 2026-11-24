from mpmath import mp, mpf, quad, sqrt, pi, log, zeta, catalan, identify, pslq
mp.dps = 90
def P(p):
    e = 1/mpf(p)
    f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
    return 8*quad(f, [0, mpf(1)/8, mpf(1)/4, mpf(3)/8, mpf(1)/2])
P0 = P(1)
def d2(h): return (P(1+h)-2*P0+P(1-h))/h**2
def d3(h): return (P(1+2*h)-2*P(1+h)+2*P(1-h)-P(1-2*h))/(2*h**3)
h1, h2 = mpf(10)**-6, mpf(10)**-6/2
D2 = (4*d2(h2)-d2(h1))/3
print("P''(1) =", mp.nstr(D2, 30)); print("(sqrt2/3) pi^2 =", mp.nstr(sqrt(2)/3*pi**2, 30)); print("diff =", mp.nstr(D2-sqrt(2)/3*pi**2, 3))
g1, g2 = mpf(10)**-5, mpf(10)**-5/2
D3 = (4*d3(g2)-d3(g1))/3
print("P'''(1) =", mp.nstr(D3, 25))
L=log(1+sqrt(2)); lam=log(2); s2=sqrt(2)
basis=[D3, pi**2, s2*pi**2, pi**2*lam, s2*pi**2*lam, pi**2*L, s2*pi**2*L, zeta(3), s2*zeta(3)]
print("pslq P''':", pslq([mp.mpf(x) for x in basis], maxcoeff=10**5, maxsteps=10**6, tol=mpf(10)**-18))
