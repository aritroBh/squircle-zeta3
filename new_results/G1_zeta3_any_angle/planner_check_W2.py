# Planner check of G1 Theorem C: W2 = -int_1^{t1} log g dt/(t+1) = zeta(2)/2 - Li2(c^2)/2, g = (t^2-1)/(2(t-cos a)), t1 = 1+2 sin(a/2)
import mpmath as mp
mp.mp.dps = 40
for a in [mp.mpf('0.7'), mp.pi/2, mp.mpf('1.3'), 2*mp.pi/3, mp.mpf('2.6')]:
    c, s, ca = mp.cos(a/2), mp.sin(a/2), mp.cos(a)
    W = -mp.quad(lambda t: mp.log((t*t-1)/(2*(t-ca)))/(t+1), [1, 1+2*s])
    print(mp.nstr(a, 6), mp.nstr(W, 25), "diff", mp.nstr(W - (mp.zeta(2)/2 - mp.polylog(2, c*c)/2), 3))
# Phi(c) = sum c^{2n}/(4n^2-1) = 1/2 - ((1-c^2)/(2c)) artanh c
c = mp.mpf('0.6'); print("Phi check", mp.nstr(mp.nsum(lambda n: c**(2*n)/(4*n*n-1), [1, mp.inf]) - (mp.mpf(1)/2 - (1-c*c)/(2*c)*mp.atanh(c)), 3))
