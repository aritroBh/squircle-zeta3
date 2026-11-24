import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, ctx
ctx.dps = 40
from step_07_partb import D_of
t0=time.monotonic()
worstD = arb(0); wc=None; n=0
x = 400.0
while x < 100000.0:
    x2 = min(100000.0, x*1.02)
    Ni = arb((x+x2)/2, (x2-x)/2)
    d = D_of(Ni)
    if d.upper() > worstD: worstD = d.upper(); wc=(x,x2)
    n+=1
    x = x2
print(f"part (b) interval cells: {n} cells, max D upper = {worstD.str(10)} at {wc}  (bound C_b=29.5)")
print("elapsed", time.monotonic()-t0)
open('progress.log','a').write(f"step: part(b) interval cells [400,1e5], max D {worstD.str(6)}\n")
