from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, im, re, sqrt, zeta,
                    gammainc, hyper, quad, log, nsum, inf, nstr, fac, si, ci, euler, psi)
mp.dps = 120
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def si2_series(z):
    return nsum(lambda k: (-1)**k * z**(2*k+1)/((2*k+1)**2 * fac(2*k+1)), [0, inf], tol=mpf('1e-110'))

# constant C = Si2(z) - (pi/2) log z limit
C1 = si2_series(mpf(1))   # ∫_0^1 Si/u
H = lambda m: quad(lambda u: (si(u)-pi/2)/u, [m*pi, (m+1)*pi])
C2 = quad(lambda u: (si(u)-pi/2)/u, [1, pi]) + nsum(lambda m: H(m), [1, inf])
C = C1 + C2
print("C =", N(C,60))
# cross-check C at z=25: si2_series(25) - (pi/2)log25 ≈ C + O(1e-3)
print("check z=25:", N(si2_series(mpf(25)) - (pi/2)*log(25) - C, 10))

# epsilon(n) for large n via asymptotic: ε(z)=Σ_{k}(-1)^{k+1}(2k)! S_{2k+3}(z), S_m = (−1)^n Σ_{j} (rising)/z^{m+2j}
def eps_asym(n, K=12, J=14):
    z = n*pi; sgn = (-1)**n
    def S(m):
        t = mpf(0); r = mpf(1)
        for j in range(J):
            term = r / z**(m+2*j)
            t += term
            r *= (m+2*j)*(m+2*j+1)
        return sgn*t
    tot = mpf(0)
    for k in range(K):
        tot += (-1)**(k+1) * fac(2*k) * S(2*k+3)
    return tot

NHEAD = 40
head = nsum(lambda n: si2_series(n*pi)/n**3, [1, NHEAD])
mid  = (pi/2)*( log(pi)*(zeta(3) - nsum(lambda n:1/n**3,[1,NHEAD])) - (zeta(3, derivative=1) + nsum(lambda n: log(n)/n**3,[1,NHEAD])) ) \
       + C*(zeta(3) - nsum(lambda n:1/n**3,[1,NHEAD]))
# note: -sum_{n>N} log n/n^3 = zeta'(3)+sum_{n<=N}log n/n^3 ; sign handled:
mid  = (pi/2)*( log(pi)*(zeta(3) - nsum(lambda n:1/n**3,[1,NHEAD])) + (zeta(3, derivative=1) + nsum(lambda n: log(n)/n**3,[1,NHEAD])) ) \
       + C*(zeta(3) - nsum(lambda n:1/n**3,[1,NHEAD]))
eps_tail = nsum(lambda n: eps_asym(n)/n**3, [NHEAD+1, inf])
S_Si2 = head + mid + eps_tail
print("S_Si2 numeric:", N(S_Si2,60))
print("49 pi^3/432  :", N(49*pi**3/432,60))
print("diff:", nstr(abs(S_Si2 - 49*pi**3/432),3))
logline("step_09: C and S_Si2 computed; compare 49 pi^3/432")
