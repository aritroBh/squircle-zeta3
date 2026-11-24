import subprocess
PY = '<project>/.venv/bin/python'
test = r'''
import mpmath as mp
mp.mp.dps = 50
# regular hexagon, K_p = {sum_j w |<n_j,x>|^p <= 1}, w = 1/sqrt(3) (edge length 1 at p=1)
# perimeter via Cauchy: P(p) = int_0^{2pi} h_p(theta) dtheta, h_p = support function
sq3 = mp.sqrt(3)
w = 1/sq3
normals = [(mp.cos(j*mp.pi/3), mp.sin(j*mp.pi/3)) for j in range(3)]
def Kp_f(x, y, p):
    return sum(abs(n[0]*x+n[1]*y)**p for n in normals)*w**p if False else sum(w**p*abs(n[0]*x+n[1]*y)**p for n in normals)
# support function: h_p(u) = max over K_p of <x,u>. Compute by solving linear program:
# K_p is convex; maximize <x,u> s.t. f(x) <= 1. Use mpmath findroot on the boundary + grid search.
def h_p(ux, uy, p):
    best = mp.mpf(0)
    # sample boundary points: f(x,y) = 1. Parametrize by angle, solve radial.
    import math
    for i in range(720):
        th = 2*mp.pi*i/720
        dx, dy = mp.cos(th), mp.sin(th)
        # find max r with f(r*dx, r*dy) <= 1
        g = lambda r: sum(w**p*abs(n[0]*dx+n[1]*dy)**p for n in normals)*r**p - 1
        # r = (sum(...))^(-1/p)
        s = sum(w**p*abs(n[0]*dx+n[1]*dy)**p for n in normals)
        r = s**(-1/mp.mpf(p))
        best = max(best, r*(dx*ux+dy*uy))
    return best
def P(p):
    return mp.quad(lambda th: h_p(mp.cos(th), mp.sin(th), p), [0, 2*mp.pi])
for p in [1, 2]:
    print("P(%d) =" % p, mp.nstr(P(p), 20))
# P(1) should be 6 (perimeter of regular hexagon with edge 1)
'''
r = subprocess.run([PY,'-c',test], capture_output=True, text=True, timeout=1200)
print(r.stdout); print(r.stderr[-600:])
