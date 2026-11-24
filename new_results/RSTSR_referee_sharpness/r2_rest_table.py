"""Part A: k = 100..400 step 25.  rest = M_k(data ball, D13f) - M_main (own Arb integral);
S_k (own Arb integral with the closed form A2); U = rest - S.  Compare with
leading law (sqrt pi/N)e^{-th}cos(th-pi/4), the varpi2 law, the C_U bound, and the referee's own U_k scale."""
from rcommon import *
setup(256)
pi, om, a = consts()
M = load_M()
sp = pi.sqrt()
v2 = acb(arb(2).sqrt()*(12+33*pi+2*pi**2), arb(2).sqrt()*(12-33*pi+2*pi**2))/(48*sp)
CU = arb("10583366.83")
print("k  rest  rad(rest)  rest/S-1  nr  cos(th-pi/4)  nr(varpi2)  |U|  |U|/CUbound  u=|U|N e^{sqrt(3piN)}/sqrt(pi)")
for k in range(100, 401, 25):
    N = k+1; th = (2*pi*N).sqrt()
    mm = 2/pi*Zmain(k).imag
    S = 2/pi*Z2(k).imag
    rest = M[k] - mm; U = rest - S
    nr = N*th.exp()*rest/sp
    pred2 = (acb(0, 1)*(th - pi/4)).exp()*(1 + v2/arb(N).sqrt())
    bound = CU*(-arb(3)/2*(pi*(N-3)).sqrt()).exp()
    u = abs(U)*N*(3*pi*N).sqrt().exp()/sp
    print(k, rest.mid().str(16), rest.rad().str(3), (rest/S-1).str(4), nr.str(8), (th-pi/4).cos().str(7),
          pred2.real.str(7), abs(U).str(4), (abs(U)/bound).str(3), u.str(4), flush=True)
