# ======== STEP 1 independent check (final): Cauchy a_n via accurate u-variable P quadrature ========
from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi, cos, sin
import json
mp.dps = 100
def P_complex(p):
    e = 1/p
    g = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
    return 8*quad(g, [log(2), 4, 10, 30, mp.inf], method='tanh-sinh')
r = mpf('0.12'); M = 96
vals = []
for k in range(M):
    vals.append(P_complex(1 + r*mpc(cos(2*pi*k/M), sin(2*pi*k/M)))/sqrt(2))
a_c = []
for n in range(0, 10):
    s = sum(vals[k]*mpc(cos(-2*pi*k*n/M), sin(-2*pi*k*n/M)) for k in range(M))/(M*r**n)
    a_c.append(s.real)
num = json.load(open('step2_num.json'))
print("n : a_n (Cauchy from true perimeter quadrature) : rel diff vs exact closed forms")
for n in range(0, 10):
    an = mpf(num['a'][n]); rd = abs(a_c[n]-an)/max(1, abs(an))
    print(f"{n}  {mp.nstr(a_c[n], 45)}   {mp.nstr(rd, 3)}")
