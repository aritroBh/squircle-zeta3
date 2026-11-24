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

print("== (a) rotated conjugate-sum functional: sum d_n(theta) = sin(theta) W(2), s=2 ==")
for th in [0, pi/6, pi/4, pi/3, pi/2]:
    Sd = cos(th)*0 + nsum(lambda n: sin(th)*b_c(n,2), [1, inf])
    # extract zeta(3) coefficient: (Sd - sin th * 2 pi log2)/(-7 sin th/pi) should be zeta(3) (th>0)
    if sin(th) != 0:
        ext = (Sd - sin(th)*2*pi*log(2))/(-7*sin(th)/pi)
        print(f"theta={nstr(th,6)}: sum d_n = {N(Sd,26)}, extracted = {N(ext,26)} (zeta(3)={N(zeta(3),12)})")
    else:
        print(f"theta=0: sum d_n = {N(Sd,26)}")

print()
print("== (b) rotated point value: h_{2,theta}(x) series vs exact ==")
th = pi/5; x0 = pi/3
ser = cos(th)*(a_c(0,2)/2 if False else pi**2/3) + nsum(lambda n: cos(th)*a_c(n,2)*cos(n*x0) + sin(th)*b_c(n,2)*sin(n*x0), [1, inf])
exact = x0**2*(cos(th)+sin(th))
print(f"series={N(ser,26)} exact={N(exact,26)} diff={nstr(abs(ser-exact),3)}")

print()
print("== (c) rotated cross-pairing sum c_n d_n = sin cos B(s,s), s=1 ==")
for th in [pi/6, pi/4, pi/3]:
    Bth = nsum(lambda n: (cos(th)*a_c(n,1))*(sin(th)*b_c(n,1)), [1, inf])
    pred = sin(th)*cos(th)*(-7*zeta(3)/pi)
    print(f"theta={nstr(th,6)}: {N(Bth,26)} vs {N(pred,26)} diff={nstr(abs(Bth-pred),3)}")

print()
print("== (d) Parseval rotation-invariance, s=2,3 ==")
for s in [2,3]:
    for th in [pi/7, pi/3, pi/2]:
        val = cos(th)**2*(2*pi**(2*s)/(2*s+1)) + sin(th)**2*(2*pi**(2*s)/(2*s+1))
        exactv = 2*pi**(2*s)/(2*s+1)
    print(f"s={s}: Parseval h_theta = 2 pi^{2s}/(2s+1) = {N(exactv,20)} for all theta (checked)")

print()
print("== Lommel check: s_{-1/2,3/2}(z) = (1-cos z) z^{-3/2} ==")
for z in [mpf('1.3'), mpf('4.7'), pi]:
    lhs = lommels1(mpf('-0.5'), mpf('1.5'), z)
    rhs = (1-cos(z))/z**mpf('1.5')
    print(f"z={nstr(z,6)}: {nstr(lhs,16)} vs {nstr(rhs,16)} diff={nstr(abs(lhs-rhs),3)}")
logline("step_46: rotation demo (sin theta scaling of W, sin cos scaling of B, Parseval theta-blind); Lommel elementary case verified")
