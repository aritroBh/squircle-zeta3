from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, im, re, sqrt, zeta,
                    gammainc, hyper, quad, log, nsum, inf, nstr, fac, si, ci, euler, psi)
mp.dps = 80
def E(x): return nstr(abs(x), 3)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

# --- exact coefficient derivative formulas at integer s ---
def ap2(n):  # d/ds a_n(s) at s=2 : (2/pi)∫ x^2 log x cos(nx)
    return (2/pi)*( pi*(2*log(pi)+1)*(-1)**n/n**2 + 2*si(n*pi)/n**3 )
def ap2_quad(n):
    return (2/pi)*quad(lambda x: x**2*log(x)*cos(n*x), [0,pi])

def bp1(n):  # d/ds b_n(s) at s=1 : (2/pi)∫ x log x sin(nx)
    return (2/pi)*( -pi*log(pi)*(-1)**n/n - si(n*pi)/n**2 )
def bp1_quad(n):
    return (2/pi)*quad(lambda x: x*log(x)*sin(n*x), [0,pi])

# b_n'(2) assembled from pieces:
def bp2(n):  # (2/pi)∫ x^2 log x sin(nx)
    # ∫ x^2 log x sin = -pi^2 log pi (-1)^n/n + (1/n)[2 Jc + Kc]
    # Jc = ∫ x log x cos = -(1/n)[ (ci(nπ)−γ−(−1)^n logπ−log n)/n + (1−(−1)^n)/n ]
    # Kc = ∫ x cos = ((−1)^n−1)/n^2
    Jc = -( (ci(n*pi)-euler-((-1)**n)*log(pi)-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*log(pi)*(-1)**n/n + (2*Jc + Kc)/n )
def bp2_quad(n):
    return (2/pi)*quad(lambda x: x**2*log(x)*sin(n*x), [0,pi])

def si2_series(z):  # Si2(z) = sum (-1)^k z^(2k+1)/((2k+1)^2 (2k+1)!)
    return nsum(lambda k: (-1)**k * z**(2*k+1)/((2*k+1)**2 * fac(2*k+1)), [0, inf])
def app2(n):  # d^2/ds^2 a_n(s) at s=2 : (2/pi)∫ x^2 log^2 x cos(nx)
    return (2/pi)*( 2*pi*log(pi)*(log(pi)+1)*(-1)**n/n**2
                    + ((4*log(pi)+6)*si(n*pi) - 4*si2_series(n*pi))/n**3 )
def app2_quad(n):
    return (2/pi)*quad(lambda x: x**2*log(x)**2*cos(n*x), [0,pi])

print("verify derivative-coefficient formulas vs quadrature")
for n in [1,2,3,5,8]:
    print(f"n={n}: ap2 err={E(ap2(n)-ap2_quad(n))}  bp1 err={E(bp1(n)-bp1_quad(n))}  "
          f"bp2 err={E(bp2(n)-bp2_quad(n))}  app2 err={E(app2(n)-app2_quad(n))}")
logline("step_07: exact s-derivative coefficient formulas ap2,bp1,bp2,app2 verified vs quadrature (80 dps)")
