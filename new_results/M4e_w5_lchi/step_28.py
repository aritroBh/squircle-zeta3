import mpmath as mp
mp.mp.dps = 80
# chi8: +1 for n = ±1 mod 8, -1 for n = ±3 mod 8, 0 even
def chi8(n):
    r = n % 8
    return {1:1, 7:1, 3:-1, 5:-1}.get(r, 0)
# 1) direct series
s = mp.nsum(lambda n: chi8(n)/mp.mpf(n)**5, [1, mp.inf])
# 2) Hurwitz zeta combination
h = mp.mpf(8)**-5 * (mp.zeta(5, mp.mpf(1)/8) + mp.zeta(5, mp.mpf(7)/8)
                     - mp.zeta(5, mp.mpf(3)/8) - mp.zeta(5, mp.mpf(5)/8))
# 3) Li5(xi8^k) form: 2*sqrt(2)*L = sum chi8(k) Li5(xi8^k)
xi = mp.exp(2*mp.pi*mp.j/8)
li = sum(chi8(k)*mp.polylog(5, xi**k) for k in range(1,8)) / (2*mp.sqrt(2))
print("series :", mp.nstr(s, 60))
print("hurwitz:", mp.nstr(h, 60))
print("li5form:", mp.nstr(li, 60))
print("series-hurwitz:", mp.nstr(abs(s-h), 8))
print("series-li5form:", mp.nstr(abs(s-li), 8))
# known check: zeta_{Q(sqrt2)}(5) = zeta(5)*L(chi8,5)? No: L-function factorization is for the Dedekind zeta:
# zeta_K(s) = zeta(s) L(chi8, s). Print both.
print("zeta(5)*L(chi8,5) =", mp.nstr(mp.zeta(5)*s, 30))
