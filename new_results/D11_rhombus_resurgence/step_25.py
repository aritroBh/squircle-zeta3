# Verify new-bp tanh-sinh across precisions at the worst point, then rerun Cauchy with it
from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi, cos, sin
import json
BP = None
def P_complex(p):
    e = 1/p
    g = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
    return 8*quad(g, [log(2), 2, 4, 6, 8, 12, 16, 24, 40, 80, mp.inf], method='tanh-sinh')
mp.dps = 100; v100 = P_complex(1 + mpc(0, mpf('0.12')))
mp.dps = 150; v150 = P_complex(1 + mpc(0, mpf('0.12')))
print("TS-new dps100 vs dps150:", mp.nstr(abs(v100-v150), 3))
mp.dps = 110
r = mpf('0.12'); M = 96
vals = [P_complex(1 + r*mpc(cos(2*pi*k/M), sin(2*pi*k/M)))/sqrt(2) for k in range(M)]
a_c = []
for n in range(0, 10):
    s = sum(vals[k]*mpc(cos(-2*pi*k*n/M), sin(-2*pi*k*n/M)) for k in range(M))/(M*r**n)
    a_c.append(s.real)
num = json.load(open('step2_num.json'))
print("\nn : rel diff Cauchy(perimeter quadrature) vs exact closed forms")
for n in range(0, 10):
    an = mpf(num['a'][n]); rd = abs(a_c[n]-an)/max(1, abs(an))
    print(f"{n}  {mp.nstr(rd, 3)}   a_n={mp.nstr(a_c[n], 35)}")
