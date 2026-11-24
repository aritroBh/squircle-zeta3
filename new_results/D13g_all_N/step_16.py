
src = open('step_38.py').read()
old = """    u = y*ureg
    logcf = ((1+u).log())/u                       # = log(1+u)/u, box-safe since u>0"""
new = """    u = y*ureg
    # logcf(u) = log(1+u)/u is DECREASING for u>0 -> tight box via point evals at endpoints
    _ulo = arb(u.lower()); _uup = arb(u.upper())
    logcf = iv((((1+_uup).log())/_uup).lower(), (((1+_ulo).log())/_ulo).upper())"""
assert old in src
src = src.replace(old, new)
open('step_38.py','w').write(src)
import importlib, sys, math, time
import step_38; importlib.reload(step_38)
from step_38 import certified_y38
ys = lambda N: math.sqrt(math.pi/(2*N))
for (n1,n2) in [(100,150),(150,225),(225,340),(340,400)]:
    r = certified_y38(ys(n2), ys(n1), nsub=4)
    print(f"[{n1},{n2}]: C={float(r['C'].upper()):.2f}  TLf*N1.5={float(r['TLf'].mid()):.3f}")
