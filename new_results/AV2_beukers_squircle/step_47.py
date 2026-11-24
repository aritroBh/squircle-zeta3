
from mpmath import mp, mpf, zeta, fabs, nstr, log, fsum
mp.dps = 260
from math import comb, factorial, log as llog
N = 600
Hs = [mpf(0)] + [fsum([mpf(1)/k for k in range(1, n+1)]) for n in range(1, 3*N+10)]
fact = [mpf(1)]
for k in range(1, 3*N+10):
    fact.append(fact[-1]*k)
S = mpf(0)
Km = []
for m in range(0, N):
    sm = mpf(0)
    Hm = Hs[m]
    for j in range(0, m+1):
        J = fact[j]**2 * fact[m]/fact[2*j+m+3]
        dJ = J*(-3 - 2*(j+1)*Hs[j] - (m+1)*Hm + (2*j+m+3)*Hs[2*j+m+3])
        sm += ((-1)**j * comb(m, j)) * dJ
    S += sm
    Km.append(sm)
print("partial sum m<%d: %s" % (N, nstr(S, 60)))
# tail model: K_m = (a ln^2 m + b ln m + c)/m^2 ; fit on m in [300, 600]
import mpmath
A = []; b = []
for m in range(300, N, 17):
    L = log(m)
    A.append([L*L, L, mpf(1)])
    b.append(Km[m]*m*m)
Am = mpmath.matrix(A); bm = mpmath.matrix(b)
sol = mpmath.lu_solve(Am.T*Am, Am.T*bm)
a, bb, c = sol[0], sol[1], sol[2]
resid = Am*sol-bm
print("tail model K_m m^2 = a ln^2 + b ln + c: a=%s b=%s c=%s resid=%s" % (
    nstr(a,10), nstr(bb,10), nstr(c,10), nstr(max(fabs(x) for x in resid),3)))
# tail sum from N to inf: int + Euler
# sum_{m>=N} (a ln^2 m + b ln m + c)/m^2
tail = a*( (log(N))**2/N + 2*log(N)/N + 2/N ) + bb*( log(N)/N + 1/N ) + c*(1/N)
# crude Euler-Maclaurin on g(m)=(a ln^2+b ln+c)/m^2: sum_{m>=N} ~ int_N^inf g + g(N)/2
def g(x):
    L = log(x)
    return (a*L*L + bb*L + c)/x**2
tail = mpmath.quad(g, [N, mpmath.inf]) + g(N)/2
print("tail estimate: %s" % nstr(tail, 20))
print("dB3/dp|_1 ~ %s" % nstr(S + tail, 40))
print("=> d1 = -dB3/dp = %s" % nstr(-(S+tail), 40))
