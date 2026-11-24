code = r'''
from squircle import *
from methodC import *
from mpmath import mp, zeta, mpf, nstr, fabs
import time
setdps(60)
def rel(a,b):
    d = fabs(a-b); return d/fabs(b) if d else mpf('1e-999')
t0=time.time()
# ZC self-consistency at p=1.001 and p=0.999
for p_ in [1.001, 0.999]:
    for s_ in [4]:
        z1 = ZC(p_, s_, K=300, R=16)
        z2 = ZC(p_, s_, K=800, R=20)
        print(f"p={p_} s={s_}: ZC(300)vs ZC(800) = {nstr(rel(z1,z2),4)}  ({time.time()-t0:.0f}s)")
        print("   ", nstr(z2, 45))
# ZA cost profile at p=1.001
t1=time.time()
za = ZA(1.001, 4)
print(f"ZA(1.001,4) = {nstr(za,45)}  ({time.time()-t1:.0f}s)")
'''
exec(code)
