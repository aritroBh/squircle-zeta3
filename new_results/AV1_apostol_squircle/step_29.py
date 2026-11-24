code = r'''
# step_25.py : tau-series sanity at p=1.25, 1.5 (uses correct c0..c6); close-out.
import mpmath as mp
mp.mp.dps = 105
def A_series(p, tol=mp.mpf('1e-100')):
    p = mp.mpf(p); s = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma(1+m/p)**2/(m**2*mp.gamma(1+2*m/p))
        s += t
        if m > 30 and abs(t) < tol: break
    return s
cs = [mp.mpf('0.548311355616075478824138388882'), mp.mpf('0.623642263728653966476642755554'),
      mp.mpf('0.335712456134004377918919527989'), mp.mpf('0.111148100693333219514134820899'),
      mp.mpf('0.0241049373555511495627640367923'), mp.mpf('0.00294602807596330098449447849793'),
      mp.mpf('-0.000258072152253727774978675619578')]
for p in ['1.25','1.5','2']:
    p_ = mp.mpf(p); tau = 1-1/p_
    approx = sum(cs[k]*tau**k for k in range(7))
    true = A_series(p_)
    print(f'p={p}: tau-series(7 terms)={mp.nstr(approx,18)} true={mp.nstr(true,18)} err={mp.nstr(abs(approx-true),3)}')
print()
print('W = c3 + 2*Sm*zeta3, Sm=(6+2pi/rt3)/9:')
W = mp.mpf('2.68290843678355888427196073142')
Sm = (6+2*mp.pi/mp.sqrt(3))/9
print('  -2*Sm*zeta3 =', mp.nstr(-2*Sm*mp.zeta(3),20), '  W+that =', mp.nstr(W-2*Sm*mp.zeta(3),20), '(=c3)')
'''
open('steps/step_25.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_25.py'],capture_output=True,text=True,timeout=1200)
print(r.stdout)
open('progress.log','a').write('step_25: tau-series sanity at p=1.25,1.5,2; close-out\n')
