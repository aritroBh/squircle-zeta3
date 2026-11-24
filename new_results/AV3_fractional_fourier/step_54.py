from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta, lommels1, hyper)
mp.dps = 80
def N(x,d=30): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
def b_c(n,s):
    if s==1: return 2*(-1)**(n+1)/n
    if s==2: return 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3)
def a_c(n,s):
    if s==1: return 2*((-1)**n-1)/(pi*n**2)
    if s==2: return 4*(-1)**n/n**2
W2 = 2*pi*log(2) - 7*zeta(3)/pi

print("== (a) rotated conjugate-sum: sum d_n(theta) = sin(theta) W(2), s=2 ==")
for th in [0, pi/6, pi/4, pi/3, pi/2]:
    Sd = nsum(lambda n: sin(th)*b_c(n,2), [1, inf])
    if sin(th) != 0:
        ext = (Sd - sin(th)*2*pi*log(2))/(-7*sin(th)/pi)
        print("theta=%s: sum d_n = %s, extracted = %s (zeta(3)=%s)" % (nstr(th,6), N(Sd,24), N(ext,24), N(zeta(3),10)))
    else:
        print("theta=0: sum d_n =", N(Sd,24))

print("== (b) rotated point value: h_{2,theta}(x) series vs exact ==")
th = pi/5; x0 = pi/3
ser = cos(th)*pi**2/3 + nsum(lambda n: cos(th)*a_c(n,2)*cos(n*x0) + sin(th)*b_c(n,2)*sin(n*x0), [1, inf])
exact = x0**2*(cos(th)+sin(th))
print("series=%s exact=%s diff=%s" % (N(ser,24), N(exact,24), nstr(abs(ser-exact),3)))

print("== (c) rotated cross-pairing sum c_n d_n = sin cos B(s,s), s=1 ==")
for th in [pi/6, pi/4, pi/3]:
    Bth = nsum(lambda n: (cos(th)*a_c(n,1))*(sin(th)*b_c(n,1)), [1, inf])
    pred = sin(th)*cos(th)*(-7*zeta(3)/pi)
    print("theta=%s: %s vs %s diff=%s" % (nstr(th,6), N(Bth,24), N(pred,24), nstr(abs(Bth-pred),3)))

print("== (d) Parseval rotation-invariance ==")
for s in [2,3]:
    print("s=%d: (1/pi) int h_theta^2 = 2 pi^%d/%d = %s for every theta" % (s, 2*s, 2*s+1, N(2*pi**(2*s)/(2*s+1),18)))

print("== Lommel check: s_{-1/2,3/2}(z) = (1-cos z) z^{-3/2} ==")
for z in [mpf('1.3'), mpf('4.7'), pi]:
    lhs = lommels1(mpf('-0.5'), mpf('1.5'), z)
    rhs = (1-cos(z))/z**mpf('1.5')
    print("z=%s: %s vs %s diff=%s" % (nstr(z,6), nstr(lhs,14), nstr(rhs,14), nstr(abs(lhs-rhs),3)))
logline("step_47: rotation demo and Lommel check done")
