import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))

# true saddle parametrized by t = rho*y:  rho^2 = (1+rt2 t+t^2)/(1+rt2 t), y^2 = t^2(1+rt2 t)/(1+rt2 t+t^2)
def saddle_from_t(t):
    rho2 = (1+sqrt2*t+t*t)/(1+sqrt2*t)
    rho = mp.sqrt(rho2); y = t/rho; N = pi/(2*y*y)
    return rho, y, N

worst = mp.mpf('9'); warg=None
worstN100 = mp.mpf('9'); wargN=None
import random
random.seed(31)
ts = [mp.mpf(10)**(k/100) for k in range(-800, 1)] + [mp.mpf(str(random.uniform(1e-8,1.0))) for _ in range(500)]
for t in ts:
    rho, y, N = saddle_from_t(t)
    x0 = 1/(rho*y); th = pi/(sqrt2*y)
    lhs = mp.re(Phi(x0,N)) + th
    mSlo = (N*t**3/2)*(sqrt2/3)/(1+sqrt2*t+t*t)
    mShi = (N*t**3/2)*(sqrt2/3)
    margin = lhs - mSlo          # must be >= 0 (and lhs <= mShi since lhs < rhs <= mShi)
    if margin < worst: worst, warg = margin, (t, N)
    if N >= 100 and margin < worstN100: worstN100, wargN = margin, (t, N)
print("mS box membership (true saddle): min(lhs - mS.lower) over t in (0,1] (must be >=0):")
print("   all t:      ", mp.nstr(worst,10), "at t,N =", [mp.nstr(v,8) for v in warg])
print("   N>=100 only:", mp.nstr(worstN100,10), "at t,N =", [mp.nstr(v,8) for v in wargN])
# also check lhs <= mShi
v2 = mp.mpf('-9')
for t in [mp.mpf(10)**(k/50) for k in range(-400, 1)]:
    rho, y, N = saddle_from_t(t)
    x0 = 1/(rho*y); th = pi/(sqrt2*y)
    lhs = mp.re(Phi(x0,N)) + th
    v2 = max(v2, lhs - (N*t**3/2)*(sqrt2/3))
print("max(lhs - mS.upper) (must be <0):", mp.nstr(v2,10))
# relative margins at the extremes used: N=100 (t~0.127), N=1e10
for Nv in [100, 400, 1e6, 1e10]:
    y = mp.sqrt(pi/(2*Nv))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2); t = rho*y; x0 = 1/(rho*y); th = pi/(sqrt2*y)
    lhs = mp.re(Phi(x0,N)) + th
    mSlo = (N*t**3/2)*(sqrt2/3)/(1+sqrt2*t+t*t); mShi = (N*t**3/2)*(sqrt2/3)
    print(f"N={Nv:g}: lhs={mp.nstr(lhs,12)}, box=[{mp.nstr(mSlo,12)},{mp.nstr(mShi,12)}], contains: {mSlo <= lhs <= mShi}")
