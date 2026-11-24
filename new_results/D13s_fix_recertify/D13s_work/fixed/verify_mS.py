
# H3 (mS): identity RePhi(x0)+th = -(N/2)S(t) - (pi rt2/(4t))(rho-1)^2  and h3-box containment
import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))
worst = mp.mpf('0'); minlo = mp.mpf('9'); minhi = mp.mpf('9'); minlo_dom = mp.mpf('9')
for k in range(-400, 401):
    t = mp.mpf(10)**(mp.mpf(k)/40)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    S = mp.log(1+sqrt2*t+t*t) - sqrt2*t
    lhs = mp.re(Phi(x0,N)) + th
    rhs = -(N/2)*S - (pi*sqrt2/(4*t))*(rho-1)**2
    worst = max(worst, abs(lhs-rhs)/abs(lhs))
    blo = (N*t**3/2)*(sqrt2/3)/(1+sqrt2*t+t*t)
    bhi = (N*t**3/2)*(sqrt2/3)
    minlo = min(minlo, lhs - blo)
    minhi = min(minhi, bhi - lhs)
    if N >= 100: minlo_dom = min(minlo_dom, lhs - blo)
print("identity worst rel resid (t in [1e-10,1e10]) =", mp.nstr(worst, 8))
print("upper box: min(box_hi - truth) over all t =", mp.nstr(minhi, 8), " (>0: valid for ALL t>0)")
print("lower box: min(truth - box_lo) over all t =", mp.nstr(minlo, 8), " (fails at large t, i.e. small N)")
print("lower box: min(truth - box_lo) for N>=100 =", mp.nstr(minlo_dom, 8), " (>0: valid on certified domain)")
