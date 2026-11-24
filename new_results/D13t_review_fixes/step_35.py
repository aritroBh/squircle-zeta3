import mpmath as mp
mp.mp.dps = 60
r2 = mp.sqrt(2)
def Hp3_mp(t):
    den = 1+r2*t+t*t
    q = mp.sqrt((1+r2*t)/den)
    w = 1/(den*(1+q))
    Gp = (r2/3)*(r2+2*t)/den**2
    hp = r2*w*w*t + r2*w*t*(2+r2*t)/(q*den**2)
    return Gp-hp
mn = mp.mpf('9'); at=None
for k in range(200001):
    t = mp.mpf(k)*mp.mpf('0.13')/200000
    v = Hp3_mp(t)
    if v < mn: mn=v; at=t
print(f"R5(a) true min of Hp3 on [0,0.13] (200001 pts) = {mp.nstr(mn,12)} at t = {mp.nstr(at,8)}")
print(f"      single-box lower 0.153836 <= true min {mp.nstr(mn,8)}: {float(mn) >= 0.1538363242986387}  (box is a valid outer enclosure, conservative)")

# (b) symbolic: H'(t) = t^3 Hp3(t) along the saddle branch?
import sympy as sp
t = sp.symbols('t', positive=True)
rt2s = sp.sqrt(2)
rho = sp.sqrt(1 + t**2/(1+rt2s*t))
psi = rt2s*t - sp.log(1+rt2s*t+t**2)
den = 1+rt2s*t+t**2
H = psi - (rt2s/3)*t**3/den - rt2s*t*(rho-1)**2/rho**2
Hp = sp.diff(H, t)
q = sp.sqrt((1+rt2s*t)/den)
w = 1/(den*(1+q))
Hp3s = (rt2s/3)*(rt2s+2*t)/den**2 - rt2s*w**2*t - rt2s*w*t*(2+rt2s*t)/(q*den**2)
diff = sp.simplify(sp.expand(Hp - t**3*Hp3s))
print("(b) sympy simplify(H' - t^3 Hp3) =", diff)
