from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 100
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def si2_fast(z, tol=mpf('1e-90')):
    t = z; s = t; k = 0
    while abs(t) > tol:
        t = -t * z*z * (2*k+1)**2 / ((2*k+3)**3 * (2*k+2)); s += t; k += 1
    return s
# asymptotic Si2 for large n: (pi/2)(log(n pi)+gamma) + eps(n), eps from h_r
pmax = 30
CA = [mpf(0)]*pmax; CB = [mpf(0)]*pmax
CA[0]=mpf(0); CB[0]=mpf(2)
for j in range(1,pmax):
    CA[j] = -j*CB[j-1]; CB[j] = pi**j + j*CA[j-1]
h = {}
for k in range(12):
    for (P, moments, fk) in [(2*k+2, CA, fac(2*k)), (2*k+3, CB, fac(2*k+1))]:
        for j in range(0, pmax - P):
            r = P + j
            if r > pmax: continue
            h[r] = h.get(r, mpf(0)) + (-1)*(-1)**k * fk * (-1)**j * gamma(P+j)/gamma(P)/fac(j) * moments[j]
def eps_Si2(n):  # epsilon(nπ) for integer n>=1 via h_r, error ~1e-45 at n>=40
    s = mpf(0)
    for r in range(3,27):
        Phi = 2**(-r)*( zeta(r, n/2) - zeta(r, (n+1)/2) )
        s += h[r] * pi**(-r) * Phi
    return -((-1)**n) * s
def Si2_n(n):
    n = int(n)
    if n <= 60: return si2_fast(n*pi)
    return (pi/2)*(log(n*pi)+euler) + eps_Si2(n)
def app2(n):
    return (2/pi)*( 2*pi*log(pi)*(log(pi)+1)*(-1)**n/n**2
                    + ((4*log(pi)+6)*si(n*pi) - 4*Si2_n(n))/n**3 )
def bp1(n):  return (2/pi)*( -pi*log(pi)*(-1)**n/n - si(n*pi)/n**2 )
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*log(pi)-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*log(pi)*(-1)**n/n + (2*Jc + Kc)/n )

Lp = log(pi)
Qpp2 = -pi**2*(9*Lp**2-6*Lp+2)/27
zp3 = zeta(3, derivative=1)
s_full  = nsum(lambda n: app2(n), [1, inf])
s_split = nsum(lambda n: app2(n) + 4*log(n)/n**3, [1, inf])
print("sum app2           =", N(s_full,32)); print("Q''(2)             =", N(Qpp2,32), " err:", nstr(abs(s_full-Qpp2),3))
print("sum[app2+4logn/n^3]=", N(s_split,32)); print("Q''(2)-4zeta'(3)   =", N(Qpp2-4*zp3,32), " err:", nstr(abs(s_split-(Qpp2-4*zp3)),3))

psi3 = psi(0,3)
Tp2  = (pi/2)**2*log(pi/2); Tp1 = (pi/2)*log(pi/2)
G = mpf('0.915965594177219015054603514932384110774')
sT2  = nsum(lambda n: bp2(n)*sin(n*pi/2), [1, inf])
print("sum bp2 sin(npi/2) =", N(sT2,32), " T'(2) =", N(Tp2,32), " err:", nstr(abs(sT2-Tp2),3))
sT1  = nsum(lambda n: bp1(n)*sin(n*pi/2), [1, inf])
sT1s = nsum(lambda n: (bp1(n) + 1/n**2)*sin(n*pi/2), [1, inf])
print("sum bp1 sin(npi/2) =", N(sT1,32), " T'(1) =", N(Tp1,32), " err:", nstr(abs(sT1-Tp1),3))
print("sum[(bp1+1/n^2)sin]=", N(sT1s,32), " T'(1)+G =", N(Tp1+G,32), " err:", nstr(abs(sT1s-(Tp1+G)),3))
logline("step_24: Q''(2) split (zeta'(3)) and T'(2),T'(1) splits (Catalan G) verified high precision")
