from mpmath import mp, mpf, sqrt, exp
mp.dps = 50
def v(u):
    return sqrt((1-exp(-2*u))/(2*u)) if u != 0 else 1
def rho3(u):
    return v(u) - 1 + u/2 - 5*u**2/24 + u**3/16
worst = 0; arg=None
for k in range(1, 2001):
    u = mpf(k)/2000
    val = abs(rho3(u))/u**3
    if val > worst: worst = val; arg = u
print(f"max |rho3(u)|/u^3 on (0,1] = {mp.nstr(worst,12)} at u={arg}")
print(f"limit u->0 = a4 = {mp.nstr(abs(rho3(mpf('0.001')))/mpf('0.001')**3, 12)}")
# also on [1, 10] for the 2.42 claim
worst2 = 0
for k in range(1000, 20001):
    u = mpf(k)/1000
    val = abs(rho3(u))/u**3
    if val > worst2: worst2 = val
print(f"max |rho3(u)|/u^3 on [1,20] = {mp.nstr(worst2,12)}  (claim <= 2.42)")
