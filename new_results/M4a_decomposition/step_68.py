import subprocess, sys
code = r'''"""Timing probe for high-precision piece quadrature at 300 dps."""
import sys, json, time
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 300
r2 = mp.sqrt(2); aN = 1+r2; bN = r2-1
omega = json.load(open("omega4.json"))
terms = {tuple(map(int, k.split(","))): [(z, mp.mpf(c[0]) + r2*mp.mpf(c[1])) for z, c in v]
         for k, v in omega.items()}
pole = {"0": mp.mpf(0), "a": aN, "-a": -aN, "b": bN, "-b": -bN}
def piece(i, j):
    ts = terms[(i, j)]
    def f(v):
        l1 = mp.log((v - bN)*(v + aN)/(4*v))
        l2 = mp.log((aN - v)*(v + bN)/(4*v))
        return l1**i * l2**j * mp.fsum(c/(v - pole[z]) for z, c in ts)
    return mp.quad(f, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
for ij in [(2,1),(3,0),(2,2),(3,1),(4,0)]:
    t0 = time.time()
    val, err = piece(*ij)
    print(f"{ij}: {time.time()-t0:.1f}s err={mp.nstr(err,3)}", flush=True)
# M4 full
import sympy as sp
from mk_integrands import v_integrand, v as vs
Rk = v_integrand(4)
fs = [(i, j, sp.lambdify(vs, R, modules="mpmath")) for (i, j), R in Rk.items()]
def fM4(x):
    l1 = mp.log((x - bN)*(x + aN)/(4*x))
    l2 = mp.log((aN - x)*(x + bN)/(4*x))
    return mp.fsum(R(x)*l1**i*l2**j for i, j, R in fs)
t0 = time.time()
M4, err = mp.quad(fM4, [bN, mp.mpf(1)/2, 1, 2, aN], error=True, maxdegree=14)
print(f"M4 full: {time.time()-t0:.1f}s err={mp.nstr(err,3)}")
'''
open("step_hp_probe.py", "w").write(code)
r = subprocess.run([sys.executable, "step_hp_probe.py"], capture_output=True, text=True, timeout=2700)
print(r.stdout)
print("ERR", r.stderr[-400:])
