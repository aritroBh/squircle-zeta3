# Check the second-order ODE formula b_n''(2) = F2_n/(1-n^2) against AV8's 50-55 digit quadrature data.
# F2 = -2 L0 s1 + 2 tan(th) s1' - L0^2 sin + 2 L0 sin  (n>=3), L0 = log|tan th|, s1 = sum beta_m sin(m th).
import json
from mpmath import mp, mpf, nsum, inf, log, pi, nstr
mp.dps = 60
def beta(m):
    if m < 0: return -beta(-m)
    if m == 1: return (log(2)-1)/2
    if m % 4 == 1: return mpf(2)/((m+1)**2*(m-1))
    return mpf(-2)/((m-1)**2*(m+1))
def c(i): return (2*i+1)*beta(2*i+1)
def F2(n):
    T1 = 2*nsum(lambda i: (beta(n-2*(2*int(i)+1))+beta(n+2*(2*int(i)+1)))/(2*i+1), [0, inf])
    j = (n-1)//2
    h = c(j) + 2*nsum(lambda r: (-1)**int(r)*c(j+int(r)), [1, inf])
    T2 = 2*h
    if n % 4 == 1:
        m = (n-1)//2; O = sum(mpf(1)/(2*i-1) for i in range(1, m//2+1)); T3 = -8*O/(n-1); T4 = mpf(4)/(n+1)
    else:
        m = (n+1)//2; O = sum(mpf(1)/(2*i-1) for i in range(1, m//2+1)); T3 = 8*O/(n+1); T4 = mpf(-4)/(n-1)
    return T1+T2+T3+T4
d = json.load(open('../AV8_squircle_sine_fourier/scratchpad/bnder60.json'))
for n in range(3, 18, 2):
    v = F2(n)/(1-n*n)
    ref = mpf(d['%d,2' % n]) if isinstance(d['%d,2' % n], str) else mpf(str(d['%d,2' % n]))
    print(n, nstr(v, 30), 'diff', nstr(v-ref, 3))
