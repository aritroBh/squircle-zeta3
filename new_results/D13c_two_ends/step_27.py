
code = r'''
from mpmath import mp, mpf, mpc, pi, sqrt, exp, log, nstr, gamma, cos, sin
import json
mp.dps = 60
R = {int(k): mpf(v) for k,v in json.load(open("mkmain_cache.json")).items()}
def load_mk(fn):
    out = {}
    for line in open(fn):
        p = line.split("\t")
        if len(p) < 2: continue
        out[int(p[0])] = mpf(p[1].strip().split(" ")[0].strip("[]"))
    return out
Mkt = load_mk("../D13b_Mk_constant/mk_deg9.tsv")
def H(s): return sqrt(pi)/4*gamma(s/2)/gamma(s/2+mpf('1.5'))
rot = exp(mpc(0,pi/4))
def RePhi(x, N):
    return -pi*x/(2*sqrt(2)) + N*log(x) - N*log(1+sqrt(2)*x+x*x)/2
def Phi(x, N):
    s = rot*x
    return mpc(0,1)*pi*s/2 + N*log(s/(1+s))
def Mkmain_quad(k, dps=55):
    old = mp.dps; mp.dps = dps
    N = mpf(k+1)
    def ig(x):
        s = rot*x
        return rot*exp(mpc(0,1)*pi*s/2)*(s/(1+s))**N*H(s)/s
    try: v = mp.quad(ig, [0, mpf('0.5'), 2, 10, 50, mp.inf])
    finally: mp.dps = old
    return (2/pi)*v.imag

# saddle approximation of I with exact x* = sqrt(2N/pi), including Gaussian and prefactor
def Mkmain_saddle(k, order=0):
    N = mpf(k+1)
    xs = sqrt(2*N/pi)
    s = rot*xs
    ph = Phi(xs, N)
    # phi''(x*) exact: d^2/dx^2 [i pi s/2 + N log(s/(1+s))] = -i N (1+2s)/(s^2(1+s)^2) * (ds/dx)^2 ; ds/dx = rot
    ph2 = -mpc(0,1)*N*(1+2*s)/(s*s*(1+s)**2)*rot**2
    ph1 = mpc(0,1)*pi*rot/2 + N*rot/(s*(1+s))  # phi'(x*) (small)
    g = sqrt(2*pi/(-ph2))
    psi = rot*H(s)/s
    I = exp(ph)*g*psi*(1 + ph1*ph1/(2*(-ph2)))   # include e^{ph1 d} linear term correction: <d^2>/2 * ph1^2
    return (2/pi)*I.imag

print("validation: M_k^main (quad) vs saddle approx vs leading law, ratio:")
for k in [50, 100, 200, 400]:
    N = mpf(k+1); th = sqrt(pi*N)
    mq = Mkmain_quad(k)
    ms = Mkmain_saddle(k)
    lead = -(sqrt(pi)/N)*exp(-th)*cos(th)
    print("  k=%3d: quad=%s saddle/quad-1=%s lead/quad-1=%s" % (
        k, nstr(mq,7), nstr(ms/mq-1,3), nstr(lead/mq-1,3)))

# RePhi concavity: RePhi''(x) < 0 for all x>0 (proved symbolically in text). Verify numerically:
def RePhi2(x,N):
    return -N/x**2 + N*x*(x+sqrt(2))/(1+sqrt(2)*x+x*x)**2
bad = 0
for N in [50, 400]:
    for j in range(1, 2000):
        x = mpf(j)/10
        if RePhi2(x,N) >= 0: bad += 1
print("\nRePhi'' sign violations over x in (0,200), N=50,400:", bad)

# exact max x0 of RePhi and RePhi(x0) expansion
from mpmath import findroot
for k in [100, 400]:
    N = mpf(k+1)
    xs = sqrt(2*N/pi)
    x0 = findroot(lambda x: -pi/(2*sqrt(2)) + N/x - N*(x+1/sqrt(2))/(1+sqrt(2)*x+x*x), xs)
    print("  k=%d: x*=%s x0=%s x0-x*=%s RePhi(x0)+sqrt(pi N)=%s" % (
        k, nstr(xs,8), nstr(x0,8), nstr(x0-xs,5), nstr(RePhi(x0,N)+sqrt(pi*N),8)))
'''
open("step_13.py","w").write(code)
import subprocess, sys
r = subprocess.run([sys.executable,"step_13.py"],capture_output=True,text=True)
print(r.stdout); print("STDERR:", r.stderr[-1200:] if r.stderr else "")
