step05 = r'''import sys, time
sys.path.insert(0, '../CX3k_review_draft8')
import polygon_engine as PE
from mpmath import mp, mpf, cos, sin, pi, tanh, sinh, cosh, nstr
import step_04 as O  # my optimized module (imports ok since __main__ guard)

mp.dps = 60
vs, ws = PE.regular_vs(3)
PF = PE.PolygonFamily(vs, ws, N=6)

# 1) compare integrand series values
maxd = mpf(0)
for i in range(1, 20):
    th = pi/6 + (pi/3)*i/20
    Fo = PF.integrand_series(th)
    Fm = O.integrand_series(th)
    d = max(abs(a-b) for a,b in zip(Fo[:6], Fm[:6]))
    maxd = max(maxd, d)
print("max |F_orig - F_mine| over 19 interior pts (coeffs 0..5):", nstr(maxd, 4))

# 2) per-level diffs, ORIGINAL style (both endpoints), one wedge, kmax=3
def run(style):
    a = pi/6; b = pi/2; mid = (a+b)/2; half = (b-a)/2
    h = mpf(2); prev = None; kmax = 3; eps = mpf(10)**(-30)
    for level in range(14):
        h = h/2
        acc = [mpf(0)]*(kmax+1)
        k = 0
        while True:
            t = k*h
            arg = (pi/2)*sinh(t)
            u = tanh(arg)
            dx = half*(pi/2)*cosh(t)/(cosh(arg)**2)*h
            if k == 0:
                ths = [mid]
            else:
                ths = [mid+half*u] if style=='paired' else [mid+half*u, mid-half*u]
            contrib = mpf(0)
            for th in ths:
                F = PF.integrand_series(th)
                for j in range(kmax+1):
                    acc[j] += (2 if style=='paired' and k>0 else 1)*dx*F[j]
                contrib = max(contrib, max(abs(dx*f) for f in F[:kmax+1]))
            k += 1
            if contrib < eps and k > 8: break
            if k > 300000: break
        diff = None if prev is None else max(abs(acc[j]-prev[j]) for j in range(kmax+1))
        print("  style=%s level %2d k=%6d diff=%s" % (style, level, k, nstr(diff,4) if diff else '-'))
        prev = acc
print("ORIGINAL style:")
run('orig')
print("PAIRED style:")
run('paired')
'''
open('step_05.py','w').write(step05)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_05.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-800:])
