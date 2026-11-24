# Cauchy coefficients a_n = (1/2 pi i) oint Q(1+d) d^{-n-1} dd,  Q = P/sqrt2, r = 0.09, M = 96 trapezoids
from mpmath import mp, mpf, mpc, quad, sqrt, pi, cos, sin
import json
mp.dps = 70
def P_complex(p):
    e = 1/p
    f = lambda t: sqrt((e*t**(e-1))**2 + (e*(1-t)**(e-1))**2)
    return 8*quad(f, [0, mpf(1)/16, mpf(1)/4, mpf(1)/2])
r = mpf('0.09'); M = 96
vals = [P_complex(1 + r*mpc(cos(2*pi*k/M), sin(2*pi*k/M))) for k in range(M)]
Q = [v/sqrt(2) for v in vals]
a_cauchy = []
for n in range(0, 10):
    s = sum(Q[k]*mpc(cos(-2*pi*k*n/M), sin(-2*pi*k*n/M)) for k in range(M))/(M*r**n)
    a_cauchy.append(s.real)
num = json.load(open('step2_num.json'))
print("n   a_n (Cauchy, P quadrature)                  a_n (exact forms)                 rel diff")
for n in range(0, 10):
    an = mpf(num['a'][n])
    rd = abs(a_cauchy[n]-an)/max(1, abs(an))
    print(f"{n}  {mp.nstr(a_cauchy[n], 40)}  {mp.nstr(rd, 3)}")
