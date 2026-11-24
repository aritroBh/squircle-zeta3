code = r'''
# step_02.py : exact corner-layer representation of the deficit, and check.
# x = exp(-eps*u) substitution in  deficit = int_0^1 (1/x) log[(1-x*yp)/(1-x)] dx  gives
# deficit(eps) = eps * int_0^inf F(eps,u) du,  F = log( (1-exp(eps*Phi(u))) / (1-exp(-eps*u)) ),
# Phi(u) = -u + log(1-exp(-u)).
import mpmath as mp
mp.mp.dps = 105

def A_series(p, tol=mp.mpf('1e-100')):
    p = mp.mpf(p); s = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma(1 + m/p)**2 / (m**2 * mp.gamma(1 + 2*m/p))
        s += t
        if m > 20 and abs(t) < tol: break
    return s

def Phi(u):
    return -u + mp.log(1 - mp.exp(-u))

def F(eps, u):
    if u == 0:
        # limit: log( ( -eps*Phi - ...)/(eps u ...) ) -> log(-Phi/u) as eps->0; for eps>0 use series
        pass
    return mp.log( (1 - mp.exp(eps*Phi(u))) / (1 - mp.exp(-eps*u)) )

def deficit_int(p):
    eps = mp.mpf(1)/p
    # integrand singular-log at u=0 (integrable), exp small at infinity
    f = lambda u: F(eps, u)
    return eps * mp.quad(f, [mp.mpf('1e-12'), mp.mpf('1e-4'), mp.mpf('1e-2'), 1, 10, 40, mp.inf])

z2 = mp.zeta(2)
for p in [2, 3, 7, 25, 100]:
    d_series = z2 - A_series(p)
    d_int = deficit_int(p)
    print(f"p={p:>4}  deficit(series)={mp.nstr(d_series,20)}  deficit(corner integral)={mp.nstr(d_int,20)}  |diff|={mp.nstr(abs(d_series-d_int),3)}")
'''
open('steps/step_02.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_02.py'],capture_output=True,text=True,timeout=900)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')
open('progress.log','a').write('step_02: corner-layer integral rep of deficit derived and verified vs series\n')
