
import os, subprocess
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
code = r'''
import mpmath as mp
mp.mp.dps = 50
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x
def TMtruth(Nv):
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    y = mp.sqrt(pi/(2*N))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-45'))
    x0 = 1/(rho*y)
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N); g0 = g(x0)
    om = mp.im(Phi1(x0,N))
    c2 = abs(mp.re(p2))/2
    w = mp.sqrt(mp.log(N)+1.0)
    d = min(w/mp.sqrt(2*c2), x0/3)
    sc = N**mp.mpf('-1.5')*mp.e**(-th)
    rt = mp.sqrt(-p2/2)
    tailR = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d + 1j*om/p2))
    tailL = mp.e**(om**2/(2*p2))*mp.sqrt(-2/p2)*(sqpi/2)*mp.erfc(rt*(d - 1j*om/p2))
    return (2/pi)*abs(g0*mp.e**Phi0*(tailR+tailL))/sc
prev = None
for Nv in [1e10, 1.2e10, 1.5e10, 2e10, 3e10, 5e10, 1e11, 3e11, 1e12]:
    v = TMtruth(Nv)
    mono = "" if prev is None else ("  (decreasing)" if v < prev else "  ** INCREASING **")
    print(f"N={Nv:g}: (2/pi)TM_truth/sc = {mp.nstr(v,12)}{mono}")
    prev = v
'''
open(os.path.join(FIX,'verify_TM_mono.py'),'w').write(code)
PY='<home>)/.venv/bin/python'
r = subprocess.run([PY, 'verify_TM_mono.py'], cwd=FIX, capture_output=True, text=True, timeout=2400)
print(r.stdout); print(r.stderr[-1000:] if r.stderr else "")
