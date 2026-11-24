import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi

# (1) y0 in [y, y(1+y^2)] for all N >= 100 (largest y = 0.12531); check ratio rho=y0/y in [1, 1+y^2]
mn = mp.mpf('9'); mx = mp.mpf('-9')
for k in range(2000):
    Nv = 100 + 400*k/1999
    y = mp.sqrt(pi/(2*Nv))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    mn = min(mn, rho); mx = max(mx, rho)
print(f"rho=y0/y over N in [100,400]: min={mp.nstr(mn,10)} max={mp.nstr(mx,10)} (must be in [1, 1+y^2]=[1,1.0157])")
# monotonicity of saddle RHS t^2(1+rt2 t)/(1+rt2 t+t^2): sympy
import sympy as sp
t = sp.symbols('t', positive=True)
F = t**2*(1+sp.sqrt(2)*t)/(1+sp.sqrt(2)*t+t**2)
dF = sp.simplify(sp.diff(F,t)*(1+sp.sqrt(2)*t+t**2)**2/t)
print("(den)^2/t * F'(t) =", sp.expand(dF), " (all coeffs positive => F strictly increasing, y0 unique)")
# (2) spot-check sweep cells 5, 27, 42 by re-running certified_y38
import subprocess
r = subprocess.run(['<home>)/.venv/bin/python','-c',
"""
import sys; sys.path.insert(0,'.')
from step_38 import certified_y38
import math
for (n1,n2,expected) in [(400.0,600.0,1255.810642465483),(math.sqrt(6634399392.564008*9951599088.846012),None,None),(9951599088.846012,10000000000.0,3868.5186155400856)]:
    y1=math.sqrt(math.pi/(2*n2)); y2=math.sqrt(math.pi/(2*n1))
    coll = certified_y38(y1,y2,nsub=2)
    print(f"cell [{n1},{n2}]: C={coll['C'].upper()} expected={expected}")
"""], cwd='research/08_theory/D13h_final_cell', capture_output=True, text=True)
print(r.stdout, r.stderr[-300:] if r.stderr else "")
