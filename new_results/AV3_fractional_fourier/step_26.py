from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 100
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
def betap(s):  # beta'(s) = -sum (-1)^k log(2k+1)/(2k+1)^s
    return -nsum(lambda k: (-1)**k*log(2*k+1)/(2*k+1)**s, [0, inf])
G = beta(2)
Lp = log(pi)

print("== T'(1) ->  S_Sisin := sum Si(n pi) sin(n pi/2)/n^2  =?  (pi^2/4) log 2 ==")
gfun = lambda n: (-1)**(n+1)*(si(n*pi) - pi/2)
resid = nsum(lambda n: (si(n*pi)-pi/2)*sin(n*pi/2)/n**2, [1, inf])
S_Sisin = (pi/2)*G + resid
print("S_Sisin =", N(S_Sisin,38))
print("(pi^2/4)log2 =", N(pi**2/4*log(2),38), " diff:", nstr(abs(S_Sisin-pi**2/4*log(2)),3))
# cross-check via hump interchange with w_n=sin(n pi/2)
N0 = 60
F = lambda m: si((m+1)*pi) - si(m*pi)
def wtail(m):  # sum_{n>m} sin(n pi/2)/n^2
    return nsum(lambda j: sin((m+j)*pi/2)/(m+j)**2, [1, inf])
head = nsum(lambda n: si(n*pi)*sin(n*pi/2)/n**2, [1, N0])
tail = si(N0*pi)*wtail(N0) + nsum(lambda m: F(m)*wtail(m), [N0, inf])
print("cross-check:", N(head+tail,38), " diff:", nstr(abs(head+tail-pi**2/4*log(2)),3))
logline("step_25: S_Sisin = sum Si(n pi)sin(n pi/2)/n^2 = (pi^2/4) log 2, two methods")
