import subprocess, os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
code = r'''
import mpmath as mp
mp.mp.dps = 80
def Emp(s):
    w = 1j*s
    return abs((mp.sqrt((1-mp.e**(-2*w))/(2*w)) - 1 + w/2)/(w*w))
vals = [(Emp(i*mp.pi/2/4000), i*mp.pi/2/4000) for i in range(4001)]
v, s = max(vals)
print("mp max |E(is)| on [0,pi/2]:", mp.nstr(v, 30), "at s =", mp.nstr(s, 10))
print("5/24 =", mp.nstr(mp.mpf(5)/24, 30))
# also on [0, 2.9] (larger slit) for safety
vals2 = [(Emp(i*mp.mpf('2.9')/4000), i*mp.mpf('2.9')/4000) for i in range(4001)]
v2, s2 = max(vals2)
print("mp max |E(is)| on [0,2.9]:", mp.nstr(v2, 20), "at s =", mp.nstr(s2, 10))
'''
with open(wd+"/step_40.py","w") as fp: fp.write(code)
import time; t0=time.time()
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_40.py"], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout[-1200:]); print("STDERR:", r.stderr[-300:]); print("elapsed", time.time()-t0)
