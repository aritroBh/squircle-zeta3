
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi
mp.dps = 60
s2 = mp.sqrt(2)

def y_p(sig, p, q):
    sp_ = mp.sqrt(p)*sig
    return mp.sqrt(1+q*q*sp_**2/4) - q*sp_/2
def x_p(sig, p, q, Z=None):
    sp_ = mp.sqrt(p)*sig
    return q*y_p(sig,p,q)/sp_
def phi_p(sig, p, q, Z):
    xp = x_p(sig,p,q)
    return p*xp + sig**-2*mp.log(1+Z/xp)
def phi2der(sig, p, q, Z):
    xp = x_p(sig,p,q)
    return sig**-2*(1/xp**2 - 1/(xp+Z)**2)
def B0(s, q):
    return 1 + s/(2*q) - mp.mpf(5)/8*s**2
def B1(s, q):
    return 1 + 3*s/(2*q) + (-mp.mpf(7)/4 + 3/(4*q**2))*s**2

def factors(Z, sig, g1):
    q = mp.sqrt(Z)
    lead = 2*(s2-1)*q/sig
    D = phi_p(sig,2,q,Z) - phi_p(sig,1,q,Z)
    F1 = exp(-(D-lead))/exp(Z/2)
    w1 = x_p(sig,1,q,Z)/y_p(sig,1,q)**2
    w2 = x_p(sig,2,q,Z)/y_p(sig,2,q)**2
    Pref = sqrt(phi2der(sig,1,q,Z)/phi2der(sig,2,q,Z))*(x_p(sig,1,q,Z)/x_p(sig,2,q,Z))*(w2/w1)**mp.mpf('1.5')
    F2 = 2*Pref
    F3 = B0(s2*sig,q)/B0(sig,q)
    N = sig**-2
    m1_1 = mp.mpf('1.5')*g1*w1/N * B1(sig,q)/B0(sig,q)
    m1_2 = mp.mpf('1.5')*g1*w2/N * B1(s2*sig,q)/B0(s2*sig,q)
    F4 = (1+m1_2)/(1+m1_1)
    return F1,F2,F3,F4

def sig1(f, Z, g1):
    # finite difference of (F-1)/sig at small sig
    s1, s2v = mp.mpf('0.002'), mp.mpf('0.001')
    v1 = (f(Z,s1,g1)-1)/s1
    v2 = (f(Z,s2v,g1)-1)/s2v
    return v2 + (v2-v1)   # richardson

for Z, g1 in [(1j*mp.pi/2, -mp.mpf('0.5')), (1j*5*mp.pi/6, None)]:
    if g1 is None:
        a = 5*mp.pi/6
        g1 = -mp.mpf(3)/4 - 1j*mp.e**(1j*a)/(4*mp.sin(a))
    q = mp.sqrt(Z)
    f1 = sig1(lambda Zz,s,g: factors(Zz,s,g)[0], Z, g1)
    f2 = sig1(lambda Zz,s,g: factors(Zz,s,g)[1], Z, g1)
    f3 = sig1(lambda Zz,s,g: factors(Zz,s,g)[2], Z, g1)
    f4 = sig1(lambda Zz,s,g: factors(Zz,s,g)[3], Z, g1)
    print('Z=%s:' % mp.nstr(Z,4))
    print('  f1(S_exp) =', mp.nstr(f1,14), '  expect -(2√2-1)q^3/12 =', mp.nstr(-(2*s2-1)*q**3/12,14))
    print('  f2(Pref)  =', mp.nstr(f2,14), '  expect 5(√2-1)q/4 =', mp.nstr(5*(s2-1)*q/4,14))
    print('  f3(Br)    =', mp.nstr(f3,14), '  expect (√2-1)/(2q) =', mp.nstr((s2-1)/(2*q),14))
    print('  f4(Mr)    =', mp.nstr(f4,14), '  expect (3/2)g1 q (1/√2-1) =', mp.nstr(mp.mpf('1.5')*g1*q*(1/s2-1),14))
    d1 = f1+f2+f3+f4
    print('  d1 total  =', mp.nstr(d1,14))
    myform = (-mp.mpf('1.5')*(1-1/s2)*g1*q + (s2-1)/(2*q) - (2*s2-1)*q**3/12 + 5*(s2-1)*q/4)
    print('  my sigma^1 formula =', mp.nstr(myform,14))
