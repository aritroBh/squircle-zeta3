import mpmath as mp
mp.mp.dps = 60
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
s = mp.mpf('0.99')
# F1 = int_0^1 x^{-s} h(x) dx. Brute force: sum over geometric slices x in [2^{-n-1}, 2^{-n}] + [1/2, 1]
tot = mp.mpf(0)
x = mp.mpf(1)
for n in range(400):
    seg = mp.quad(lambda t: t**(-s)*h(t), [x/2, x])
    tot += seg
    x /= 2
    if x < mp.mpf(10)**(-50): break
F1 = tot
print("F1 slices:", mp.nstr(F1,15))
# int_0^1 x^-s dx = 100 check of slice scheme
tot2 = mp.mpf(0); x = mp.mpf(1)
for n in range(400):
    tot2 += mp.quad(lambda t: t**(-s), [x/2, x])
    x /= 2
    if x < mp.mpf(10)**(-50): break
print("int x^-s via slices:", mp.nstr(tot2,15), "vs 100 err", mp.nstr(abs(tot2-100),3))
# F1 by series: h(x) = sum a_n x^n? h = 1 - x/2 - x^2/2 + 3x^3/8 + 3x^4/8 - 5x^5/16 - 5x^6/16 + 35x^7/128...
# general: h(x) = sum_{n>=1} (-1)^n c_n x^{2n-1} + sum_{n>=0} (-1)^n c_n x^{2n}, c_n = binom(2n,n)/4^n
cn = lambda n: mp.binomial(2*n,n)/4**n
F1s = mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n+1-s), [0, mp.inf]) + mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n-s), [1, mp.inf])
print("F1 partial-fraction series:", mp.nstr(F1s,15))
F2s = mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n+1+s), [0, mp.inf]) + mp.nsum(lambda n: (-1)**int(n)*cn(n)/(2*n+s), [1, mp.inf])
print("F1+F2 series:", mp.nstr(F1s+F2s,15))
rhs = (mp.gamma(1+s/2)*mp.gamma(mp.mpf('0.5')-s/2) - mp.gamma(1-s/2)*mp.gamma(mp.mpf('0.5')+s/2))/(mp.sqrt(mp.pi)*s)
print("rhs:", mp.nstr(rhs,15), " series-rhs:", mp.nstr(abs(F1s+F2s-rhs),3))
