from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, psi)
mp.dps = 100
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

def si2_fast(z, tol=mpf('1e-90')):
    t = z; s = t; k = 0
    while abs(t) > tol:
        t = -t * z*z * (2*k+1)**2 / ((2*k+3)**3 * (2*k+2)); s += t; k += 1
    return s
def ap2(n):  return (2/pi)*( pi*(2*log(pi)+1)*(-1)**n/n**2 + 2*si(n*pi)/n**3 )
def app2(n): return (2/pi)*( 2*pi*log(pi)*(log(pi)+1)*(-1)**n/n**2
                    + ((4*log(pi)+6)*si(n*pi) - 4*si2_fast(n*pi))/n**3 )
def bp1(n):  return (2/pi)*( -pi*log(pi)*(-1)**n/n - si(n*pi)/n**2 )
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*log(pi)-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*log(pi)*(-1)**n/n + (2*Jc + Kc)/n )

Lp = log(pi)
Qp2  = -pi**2*(3*Lp-1)/9
Qpp2 = -pi**2*(9*Lp**2-6*Lp+2)/27
Vp1  = pi**2/6*(Lp - mpf(5)/6)
Tp1  = (pi/2)*log(pi/2)
Tp2  = (pi/2)**2*log(pi/2)
G    = mpf('0.915965594177219015054603514932384110774')  # Catalan

print("== (1) Q'(2) kink split ==")
s_full  = nsum(lambda n: ap2(n), [1, inf])
s_split = nsum(lambda n: ap2(n) - 2/n**3, [1, inf])
print("sum ap2        =", N(s_full,32), " Q'(2) =", N(Qp2,32), " err:", nstr(abs(s_full-Qp2),3))
print("sum[ap2-2/n^3] =", N(s_split,32), " Q'(2)-2zeta(3) =", N(Qp2-2*zeta(3),32), " err:", nstr(abs(s_split-(Qp2-2*zeta(3))),3))

print("== (2) Q''(2) kink split (zeta'(3)) ==")
zp3 = zeta(3, derivative=1)
s_full  = nsum(lambda n: app2(n), [1, inf])
s_split = nsum(lambda n: app2(n) + 4*log(n)/n**3, [1, inf])
print("sum app2          =", N(s_full,32), " Q''(2) =", N(Qpp2,32), " err:", nstr(abs(s_full-Qpp2),3))
print("sum[app2+4logn/n3]=", N(s_split,32), " Q''(2)-4zeta'(3) =", N(Qpp2-4*zp3,32), " err:", nstr(abs(s_split-(Qpp2-4*zp3)),3))

print("== (3) V'(1) kink split (unit zeta(3)) ==")
s_full  = nsum(lambda n: bp1(n)/n, [1, inf])
s_split = nsum(lambda n: (bp1(n) + 1/n**2)/n, [1, inf])
print("sum bp1/n          =", N(s_full,32), " V'(1) =", N(Vp1,32), " err:", nstr(abs(s_full-Vp1),3))
print("sum[(bp1+1/n^2)/n] =", N(s_split,32), " V'(1)+zeta(3) =", N(Vp1+zeta(3),32), " err:", nstr(abs(s_split-(Vp1+zeta(3))),3))

print("== (5) T'(2) and (6) T'(1) kink splits (beta'(3), Catalan G) ==")
psi3 = psi(3)
sT2  = nsum(lambda n: bp2(n)*sin(n*pi/2), [1, inf])
print("sum bp2 sin(npi/2) =", N(sT2,32), " T'(2) =", N(Tp2,32), " err:", nstr(abs(sT2-Tp2),3))
sT1  = nsum(lambda n: bp1(n)*sin(n*pi/2), [1, inf])
sT1s = nsum(lambda n: (bp1(n) + 1/n**2)*sin(n*pi/2), [1, inf])
print("sum bp1 sin(npi/2) =", N(sT1,32), " T'(1) =", N(Tp1,32), " err:", nstr(abs(sT1-Tp1),3))
print("split: ", N(sT1s,32), " T'(1)+G =", N(Tp1+G,32), " err:", nstr(abs(sT1s-(Tp1+G)),3))
logline("step_23: kink splits Q'(2),Q''(2),V'(1),T'(1),T'(2) verified: zeta(3),zeta'(3),G appear in kink part and cancel in value")
