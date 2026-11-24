from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, zeta, si, ci, euler, psi, lommels1)
mp.dps = 100
def N(x,d=30): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
Lp = log(pi)
print("== Lommel: s_{-1/2,3/2}(z) = (1 - cos z - z sin z) z^{-3/2} ==")
for z in [mpf('1.3'), mpf('4.7'), mpf('9.2')]:
    lhs = lommels1(mpf('-0.5'), mpf('1.5'), z)
    rhs = (1-cos(z)-z*sin(z))/z**mpf('1.5')
    print("z=%s diff=%s" % (nstr(z,6), nstr(abs(lhs-rhs),3)))

print("== T'(1) kink split (Catalan G), reindexed ==")
def bp1(n):  return (2/pi)*( -pi*Lp*(-1)**n/n - si(n*pi)/n**2 )
beta = lambda s: (zeta(s, mpf(1)/4) - zeta(s, mpf(3)/4))/4**s
G = beta(2); Tp1 = (pi/2)*log(pi/2)
sT1s = nsum(lambda k: (bp1(2*k+1) + 1/(2*k+1)**2)*(-1)**k, [0, inf])
print("sum[(bp1+1/n^2) sin] =", N(sT1s,36))
print("T'(1) + G              =", N(Tp1+G,36), " diff:", nstr(abs(sT1s-(Tp1+G)),3))

print("== W'(2) kink split ==")
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
psi3 = psi(0,3)
zp3 = zeta(3, derivative=1)
CiS = nsum(lambda n: ci(n*pi)/n**3, [1, inf])
Wp2 = 2*pi*Lp*log(2) - (4/pi)*(CiS + zp3 + Lp*3*zeta(3)/4 - euler*zeta(3)) - (21/(2*pi))*zeta(3)
def kinkp(n): return (4/pi)*(log(n)-psi3)/n**3
lhs = nsum(lambda n: bp2(n)-kinkp(n), [1, inf])
rhs = Wp2 + (4/pi)*(zp3 + psi3*zeta(3))
print("sum[bp2 - kink'] =", N(lhs,32))
print("W'(2)+(4/pi)(zeta'(3)+psi(3)zeta(3)) =", N(rhs,32), " diff:", nstr(abs(lhs-rhs),3))

print("== (c) rotated cross-pairing via partial sums ==")
def a_c1(n): return 2*((-1)**n-1)/(pi*n**2)
def b_c1(n): return 2*(-1)**(n+1)/n
M = 200000
for th in [pi/6, pi/4, pi/3]:
    Bth = sin(th)*cos(th)*nsum(lambda n: a_c1(n)*b_c1(n), [1, M])
    pred = sin(th)*cos(th)*(-7*zeta(3)/pi)
    print("theta=%s: %s vs %s diff=%s" % (nstr(th,6), N(Bth,22), N(pred,22), nstr(abs(Bth-pred),3)))
logline("step_48: Lommel fixed; T'(1) G-split; W'(2) kink split; rotated B partial sums")
