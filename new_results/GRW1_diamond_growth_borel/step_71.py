import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
# (A) Certified E* = max_{s in [0, pi/2]} |E(is)|, E(w) = [sqrt((1-e^{-2w})/(2w)) - 1 + w/2]/w^2
from flint import acb, arb, ctx
import mpmath as mp
mp.mp.dps = 60
ctx.prec = 180
def E(w):
    return (((-2*w).exp().neg()+1)/(2*w)).sqrt() .sub(1).add(w/2)/(w*w)
def E2(w):
    return (((1-(-2*w).exp())/(2*w)).sqrt() - 1 + w/2)/(w*w)
# grid + crude Lipschitz via complex intervals
import math
N = 4000
mx = arb(0); argmx = None
for i in range(N+1):
    s = arb(i*math.pi/2/N)
    w = acb(0,1)*acb(s)
    v = abs(E2(w))
    if v > mx: mx = v; argmx = float(s)
print("certified-ish max |E(is)| on grid (pi/2, 4001 pts):", mx, "at s=", argmx)
# interval check near max with a ball
s0 = arb(argmx)
ball = acb(0,1)*acb(s0).add_error(arb('1e-6'))
v = abs(E2(ball))
print("ball value at argmax +/-1e-6:", v)
# also global bound: |E| on [0, pi/2] continuous; sample derivative bound by secant
mp.mp.dps = 80
def Emp(s):
    w = 1j*s
    return abs((mp.sqrt((1-mp.e**(-2*w))/(2*w)) - 1 + w/2)/(w*w))
vals = [Emp(i*mp.pi/2/2000) for i in range(2001)]
print("mp max:", max(vals))
print("E(0+) should be 5/24 =", 5/24, " E(pi/2)=", Emp(mp.pi/2))
'''
with open(wd+"/step_39.py","w") as fp: fp.write(code)
t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_39.py"], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout[-2000:]); print("STDERR:", r.stderr[-500:]); print("elapsed", time.time()-t0)
