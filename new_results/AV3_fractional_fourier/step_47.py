from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, log, nsum, inf, nstr, fac, si, ci, zeta, euler, quad)
mp.dps = 100
def N(x,d=30): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
Lp = log(pi)
def si2_fast(z, tol=mpf('1e-90')):
    t = z; s = t; k = 0
    while abs(t) > tol:
        t = -t * z*z * (2*k+1)**2 / ((2*k+3)**3 * (2*k+2)); s += t; k += 1
    return s
def ap2(n):  return (2/pi)*( pi*(2*Lp+1)*(-1)**n/n**2 + 2*si(n*pi)/n**3 )
def bp2(n):
    Jc = -( (ci(n*pi)-euler-((-1)**n)*Lp-log(n))/n + (1-(-1)**n)/n )/n
    Kc = ((-1)**n - 1)/n**2
    return (2/pi)*( -pi**2*Lp*(-1)**n/n + (2*Jc + Kc)/n )
def bp3(n):  return (2/pi)*( -pi**3*Lp*(-1)**n/n + 3*pi*(2*Lp+1)*(-1)**n/n**3 + 6*si(n*pi)/n**4 + 2*pi*(-1)**n/n**3 )
def b_c(n,s):
    if s==2: return 2*pi*(-1)**(n+1)/n - 4*(1-(-1)**n)/(pi*n**3)
    if s==3: return 2*pi**2*(-1)**(n+1)/n + 12*(-1)**n/n**3
def a_c(n,s):
    if s==2: return 4*(-1)**n/n**2
    if s==3: return 6*pi*(-1)**n/n**2 + 12*(1-(-1)**n)/(pi*n**4)
    if s==4: return (-1)**n*(8*pi**2/n**2 - 48/n**4)
def ap3(n): return -(3*bp2(n) + b_c(n,2))/n
def ap4(n): return -(4*bp3(n) + b_c(n,3))/n
def bp4(n): return (2/pi)*(-pi**4*Lp*(-1)**n/n) + (4*ap3(n) + a_c(n,3))/n
# verify recursions vs quadrature
qa3 = max(abs(ap3(n) - (2/pi)*quad(lambda x: x**3*log(x)*cos(n*x), [0,pi])) for n in [1,2,5])
qa4 = max(abs(ap4(n) - (2/pi)*quad(lambda x: x**4*log(x)*cos(n*x), [0,pi])) for n in [1,2,5])
qb4 = max(abs(bp4(n) - (2/pi)*quad(lambda x: x**4*log(x)*sin(n*x), [0,pi])) for n in [1,2,5])
print("recursion errors: ap3", nstr(qa3,3), " ap4", nstr(qa4,3), " bp4", nstr(qb4,3))

checks = [
 ("Q'(3)", nsum(lambda n: ap3(n), [1, inf]), pi**3*(1-4*Lp)/16),
 ("P'(3)", nsum(lambda n: (-1)**n*ap3(n), [1, inf]), pi**3*(12*Lp+1)/16),
 ("Q'(4)", nsum(lambda n: ap4(n), [1, inf]), pi**4*(1-5*Lp)/25),
 ("V'(3)", nsum(lambda n: bp3(n)/n, [1, inf]), pi**4*(20*Lp-9)/400),
 ("V'(4)", nsum(lambda n: bp4(n)/n, [1, inf]), pi**5*(30*Lp-11)/900),
 ("Psi'(3)", 2*nsum(lambda n: b_c(n,3)*bp3(n), [1, inf]), 4*pi**6*(7*Lp-1)/49),
 ("Phi'(3)", 2*nsum(lambda n: a_c(n,3)*ap3(n), [1, inf]), 3*pi**6*(-5+84*Lp)/784),
 ("Phi'(4)", 2*nsum(lambda n: a_c(n,4)*ap4(n), [1, inf]), 16*pi**8*(-11+180*Lp)/10125),
 ("T'(4)", nsum(lambda k: bp4(2*k+1)*(-1)**k, [0, inf]), (pi/2)**4*log(pi/2)),
]
for name, got, exact in checks:
    print(f"{name}: series={N(got,28)}  exact={N(exact,28)}  diff={nstr(abs(got-exact),3)}")
logline("step_42: series-side verification of first derivatives at s=3,4 (Q',P',V',Psi',Phi',T')")
