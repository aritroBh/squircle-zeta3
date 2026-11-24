from mpmath import mp, mpf, nstr, pi, sin, cos, legenp, quad
mp.dps = 55
def Jfull_series3(s, ca, N=3000):
    x = -ca
    P0, P1 = mpf(1), x
    tot = P0*(1/s + 1/(1-s)) + P1*(1/(1+s) + 1/(2-s))
    for n in range(1, N-1):
        P0, P1 = P1, ((2*n+1)*x*P1 - n*P0)/(n+1)
        tot += P1*(1/(n+1+s) + 1/(n+2-s))
    return tot
def Jfull_trig(s, alpha):
    return quad(lambda u: sin(u)**(s-1)*sin(alpha+u)**(-s), [0, pi-alpha])
alpha = pi/3; ca = cos(alpha); s = mpf('0.5')
print("series:", nstr(Jfull_series3(s, ca), 25))
print("trig-quad:", nstr(Jfull_trig(s, alpha), 25))
print("candidate pi P_{s-1}(-ca)/sin(ps):", nstr(pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s), 25))
# also alternative candidate: pi P_{s-1}(+ca)/sin? or with different argument
print("alt +ca:", nstr(pi*legenp(s-1, 0, ca, type=2)/sin(pi*s), 25))
