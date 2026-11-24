code = r'''
# step_18.py : diagnose d9 residual
import mpmath as mp
mp.mp.dps = 105
def A_series_s(p, s=0, tol=None):
    p = mp.mpf(p); s = mp.mpf(s); tol = tol or mp.mpf('1e-100')
    tot = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma((m+s)/p)**2 / (p**2 * mp.gamma(1 + 2*(m+s)/p))
        tot += t
        if m > 30 and abs(t) < tol: break
    return tot
d1 = mp.mpf('1.73300259017088785677509047603076917317936436')
d2 = -mp.zeta(2)/2; d3 = mp.zeta(3)/6
d5 = (mp.zeta(2)*mp.zeta(3)-3*mp.zeta(5))/60
d7 = mp.zeta(7)/14 - mp.zeta(2)*mp.zeta(5)/42 - mp.zeta(2)**2*mp.zeta(3)/140
for p in [40, 50, 80, 120, 200, 400]:
    e = mp.mpf(1)/p
    D = mp.zeta(2) - A_series_s(p)
    R = (D - d1*e - d2*e**2 - d3*e**3 - d5*e**5 - d7*e**7)/e**9
    print(f'p={p}: D={mp.nstr(D,12)}, R={mp.nstr(R,12)}')
'''
open('steps/step_18.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_18.py'],capture_output=True,text=True,timeout=1200)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
