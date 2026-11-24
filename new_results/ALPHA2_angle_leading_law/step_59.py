from mpmath import mp, mpf, nstr, pi, sin, cos, legenp, nsum, inf
mp.dps = 55
def Jfull_series(s, ca):
    # Q^{-1} = sum P_n(-ca) w^n; Jfull = sum P_n(-ca)[1/(n+s)+1/(n+1-s)]
    return nsum(lambda n: legenp(n, 0, -ca, type=2)*(1/(n+s) + 1/(n+1-s)), [0, inf], tol=mpf('1e-50'))
for alpha in [pi/2, pi/3, 2*pi/3, pi/6]:
    ca = cos(alpha)
    for s in [mpf('0.9'), mpf('0.3'), mpf('0.5')]:
        a = Jfull_series(s, ca); b = pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s)
        print(f"alpha={nstr(alpha,6)}, s={nstr(s,4)}: |series-cand|={nstr(abs(a-b),3)}")
