import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
step45 = r'''
import json
from mpmath import mp, mpf, zeta, fabs, nstr, matrix, lu_solve
mp.dps = 100
exec(open('av2core.py').read())
print("=== B3 tau-grid near p=1 (polar 2D, two h levels) ===")
grid = {}
for p in [0.9, 0.95, 1.05, 1.1]:
    v1 = B3_tsB(p, '0.055', 180)
    v2 = B3_tsB(p, '0.045', 220)
    grid[p] = v1
    print("p=%5.2f B3=%s  |d|=%s" % (p, nstr(v1, 40), nstr(fabs(v1-v2), 3)))
B3_1 = mpf('0.24661730238492374056516806883959205599718610151827')
grid[mpf(1)] = B3_1
pts = [(mpf(1)-1/p, v) for p, v in sorted(grid.items())]
print(pts)
A = []; b = []
for (tau, v) in pts:
    A.append([mpf(1), tau, tau**2, tau**3])
    b.append(v)
Am = matrix(A); bm = matrix(b)
sol = lu_solve(Am.T*Am, Am.T*bm)
res = Am*sol - bm
print("residual:", nstr(max(fabs(x) for x in res),5))
for nm, s in zip(['b0','d1','d2','d3'], sol):
    print(nm, "=", nstr(s, 30))
print("=== B2 validation table: series vs quad2 at 120 dps ===")
mp.dps = 120
for p in [1, 1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 96]:
    s, nn = B2_series(p)
    q2 = B2_quad2(p)
    print("p=%5s |series-quad2|=%s (terms %d)" % (p, nstr(fabs(s-q2),3), nn))
'''
open('step_45.py','w').write(step45)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_45.py'],
                   capture_output=True, text=True, timeout=1750)
print(r.stdout); print("STDERR:", r.stderr[-500:] if r.returncode!=0 else "(ok)")
with open('progress.log','a') as f:
    f.write("step_45: B3 tau-grid near p=1 (d1 fit); B2 series-vs-quad2 residuals to 120 dps. DONE\n")
