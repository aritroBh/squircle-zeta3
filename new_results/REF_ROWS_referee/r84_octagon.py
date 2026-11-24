# R84: regular 2m-gon smoothing K_p = {sum_i |<n_i,x>|^p <= 1}, n_i=(cos(i pi/m), sin(i pi/m)).
# P(1+d) via truncated Taylor arithmetic in d under the polar perimeter integral.
import sys
from mpmath import mp, mpf, sqrt, log, cos, sin, pi, quad, zeta, polylog, nstr, exp, fsum
mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 45
N = 4  # orders d^0..d^3
def smul(a, b): return [fsum(a[i]*b[k-i] for i in range(k+1)) for k in range(N)]
def sexp(a):  # exp of series
    e0 = exp(a[0]); r = [e0] + [mpf(0)]*(N-1)
    # r' = a' r
    for k in range(1, N):
        r[k] = fsum(j*a[j]*r[k-j] for j in range(1, k+1))/k
    return r
def slog(a):
    r = [log(a[0])] + [mpf(0)]*(N-1)
    for k in range(1, N):
        r[k] = (a[k] - fsum(j*r[j]*a[k-j] for j in range(1, k))/k)/a[0]
    return r
def sinv(a):
    r = [1/a[0]] + [mpf(0)]*(N-1)
    for k in range(1, N):
        r[k] = -fsum(a[j]*r[k-j] for j in range(1, k+1))/a[0]
    return r
def ssqrt(a): return sexp([x/2 for x in slog(a)])
def integrand(th, m):
    S = [mpf(0)]*N; Sp = [mpf(0)]*N
    for i in range(1, m+1):
        ph = th - i*pi/m
        c = cos(ph); s = sin(ph); ac = abs(c); L = log(ac); sg = 1 if c > 0 else -1
        # |c|^(1+d) = ac * exp(d L)
        e = [ac*L**k/[1,1,2,6][k] for k in range(N)]
        S = [S[k]+e[k] for k in range(N)]
        # d/dth |c|^p = -p |c|^(p-1) sg s  ;  p|c|^{p-1} = (1+d) exp(d L)
        f = [L**k/[1,1,2,6][k] for k in range(N)]
        f = smul([mpf(1), mpf(1), mpf(0), mpf(0)], f)
        Sp = [Sp[k] - sg*s*f[k] for k in range(N)]
    invp = [mpf((-1)**k) for k in range(N)]       # 1/(1+d)
    logr = [-x for x in smul(invp, slog(S))]
    r = sexp(logr)
    rr = [-x for x in smul(invp, smul(Sp, sinv(S)))]  # r'/r
    one = [mpf(1)] + [mpf(0)]*(N-1)
    return smul(r, ssqrt([one[k] + smul(rr, rr)[k] for k in range(N)]))
def coeffs(m):
    # integrate over one edge (0, pi/m), multiply by 2m edges
    res = []
    cache = {}
    def F(th):
        key = th
        if key not in cache: cache[key] = integrand(th, m)
        return cache[key]
    out = []
    for k in range(N):
        out.append(2*m*quad(lambda t: F(t)[k], [0, pi/(2*m), pi/m]))
    return out  # [P(1), P'(1), P''(1)/2, P'''(1)/6]
if __name__ == "__main__":
    c2 = coeffs(2)
    print("diamond P(1)", nstr(c2[0], 30), " 4sqrt2:", nstr(4*sqrt(2), 30))
    print("diamond P'(1)", nstr(c2[1], 10))
    print("diamond P''(1)", nstr(2*c2[2], 30), " pi^2 sqrt2/3:", nstr(pi**2*sqrt(2)/3, 30))
    print("diamond P'''(1)", nstr(6*c2[3], 30), " 6sqrt2(z3-pi^2/3):", nstr(6*sqrt(2)*(zeta(3)-pi**2/3), 30))
    c4 = coeffs(4)
    s2 = sqrt(2); u = sqrt(2+s2)
    print("octagon P(1)", nstr(c4[0], 35), " printed 8(s2-1)sqrt(2-s2):", nstr(8*(s2-1)*sqrt(2-s2), 35))
    print("octagon P'(1)", nstr(c4[1], 35), " printed:", nstr(u*((40-28*s2)*log(2)+(24-16*s2)*log(1+s2)), 35))
    print("octagon P''(1)", nstr(2*c4[2], 35))
    print("octagon P'''(1)", nstr(6*c4[3], 40))
