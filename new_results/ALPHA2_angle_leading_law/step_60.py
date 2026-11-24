from mpmath import mp, mpf, nstr, pi, sin, cos, legenp
mp.dps = 55
def Jfull_series2(s, ca, N=200000):
    tot = mpf(0)
    for n in range(N):
        tot += legenp(n, 0, -ca, type=2)*(1/(n+s) + 1/(n+1-s))
    # tail ~ sum P_n (2n+1)/n^2 ~ n^{-3/2}; add rough tail estimate? just print
    return tot
for alpha, s in [(pi/2, mpf('0.5')), (pi/2, mpf('0.3'))]:
    ca = cos(alpha)
    a = Jfull_series2(s, ca)
    b = pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s)
    print(f"alpha={nstr(alpha,5)}, s={nstr(s,3)}: partial={nstr(a,20)} cand={nstr(b,20)} |diff|={nstr(abs(a-b),3)}")
