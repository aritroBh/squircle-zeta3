import mpmath as mp
mp.mp.dps = 50
m = 2; A = mp.pi/m
def S(th, d):
    p = 1+d
    return mp.fsum([mp.sin(th + i*A)**p for i in range(m)])
def rho_num(th, d):
    p = 1+d
    return S(th, d)**(-1/p)
def rho_th(th, d):
    h2 = mp.mpf('1e-12')
    return (rho_num(th-2*h2,d) - 8*rho_num(th-h2,d) + 8*rho_num(th+h2,d) - rho_num(th+2*h2,d))/(12*h2)
def Dens(th, d):
    r = rho_num(th,d); rt = rho_th(th,d)
    return mp.sqrt(r*r + rt*rt)

th = mp.mpf('0.3'); h = mp.mpf('1e-9')
g = lambda dd: Dens(th, dd)
d1 = (g(-2*h) - 8*g(-h) + 8*g(h) - g(2*h))/(12*h)
d2 = (-g(-2*h) + 16*g(-h) - 30*g(0) + 16*g(h) - g(2*h))/(12*h*h)
d3 = (g(-3*h) - 12*g(-2*h) + 39*g(-h) - 56*g(0) + 39*g(h) - 12*g(2*h) + g(3*h))/(6*h**3)
print("D:", mp.nstr(g(0),25))
print("D_d numeric:", mp.nstr(d1,20))
print("D_dd/2 numeric:", mp.nstr(d2/2,20))
print("D_ddd/6 numeric:", mp.nstr(d3/6,20))

# now my formula pieces
def mycoeffs(th):
    ls = [mp.log(mp.sin(th + i*A)) for i in range(m)]
    sv = [mp.sin(th + i*A) for i in range(m)]
    cv = [mp.cos(th + i*A) for i in range(m)]
    T = [mp.fsum([sv[i]*ls[i]**k for i in range(m)]) for k in range(4)]
    Tp = [mp.fsum(cv)] + [mp.fsum([cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,4)]
    S0 = T[0]; S0p = Tp[0]
    U1 = T[1]/S0
    U2 = T[2]/(2*S0) - T[1]**2/(2*S0**2)
    U3 = T[3]/(6*S0) - T[1]*T[2]/(2*S0**2) + T[1]**3/(3*S0**3)
    U1p = Tp[1]/S0 - T[1]*S0p/S0**2
    U2p = Tp[2]/(2*S0) - T[2]*S0p/(2*S0**2) - T[1]*Tp[1]/S0**2 + T[1]**2*S0p/S0**3
    U3p = (Tp[3]/(6*S0) - T[3]*S0p/(6*S0**2)
           - (Tp[1]*T[2]/(2*S0**2) + T[1]*Tp[2]/(2*S0**2) - T[1]*T[2]*S0p/S0**3)
           + (T[1]**2*Tp[1]/S0**3 - T[1]**3*S0p/S0**4))
    lS0 = mp.log(S0)
    R1 = lS0 - U1; R2 = -lS0 + U1 - U2; R3 = lS0 - U1 + U2 - U3
    Q0 = S0p/S0
    Q1 = Q0 - U1p; Q2 = -Q0 + U1p - U2p; Q3 = Q0 - U1p + U2p - U3p
    rho1 = R1; rho2 = R2 + R1**2/2; rho3 = R3 + R1*R2 + R1**3/6
    W1 = 2*Q0*Q1; W2 = 2*Q0*Q2 + Q1**2; W3 = 2*Q0*Q3 + 2*Q1*Q2
    B = 1 + Q0**2
    E1 = W1/(2*B); E2 = W2/(2*B) - W1**2/(8*B**2); E3 = W3/(2*B) - W1*W2/(4*B**2) + W1**3/(16*B**3)
    D0 = (1/S0)*mp.sqrt(B)
    return [D0, D0*(rho1 + E1), D0*(rho2 + rho1*E1 + E2), D0*(rho3 + rho2*E1 + rho1*E2 + E3)]
mc = mycoeffs(th)
print("my D0..D3:", [mp.nstr(v,20) for v in mc])
