import mpmath as mp
mp.mp.dps = 50

m = 2; A = mp.pi/m
def S(th, d):
    p = 1+d
    return mp.fsum([mp.sin(th + i*A)**p for i in range(m)])
def rho_num(th, d):
    p = 1+d
    return S(th, d)**(-1/p)

th = mp.mpf('0.3')
# numeric d-derivatives of rho at d=0 via 4th-order FD with tiny h
h = mp.mpf('1e-8')
f = lambda dd: rho_num(th, dd)
d1 = (f(-2*h) - 8*f(-h) + 8*f(h) - f(-2*h))/(12*h)
# analytic pieces from my formulas
sv = [mp.sin(th+i*A) for i in range(m)]
ls = [mp.log(v) for v in sv]
T = [mp.fsum([sv[i]*ls[i]**k for i in range(m)]) for k in range(4)]
S0 = T[0]
U1 = T[1]/S0
R1 = mp.log(S0) - U1
print("numeric d rho/dd at d=0:", mp.nstr(d1, 30))
print("rho0*R1               :", mp.nstr((1/S0)*R1, 30))
print("rho0                  :", mp.nstr(1/S0, 30))
print("R1 =", mp.nstr(R1, 30))
# check: log rho = -e log S, e = 1/p. d/dd log rho = -(de/dd) log S - e d/dd log S
# at d=0: de/dd = -1; d/dd log S = U1. so R1 = log S0 - U1. rho1 = rho0*R1. Verify separately:
p = 1
dlogS = U1
print("direct: -(de/dd)*logS - e*dlogS =", mp.nstr(mp.log(S0) - U1, 30))
