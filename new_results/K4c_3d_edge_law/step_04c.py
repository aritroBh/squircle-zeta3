"""step_04c.py -- cubo inner-refinement trend."""
import os, json, time
os.chdir("research/08_theory/K4c_3d_edge_law")
import mpmath as mp
from mpmath import mpf
mp.mp.dps = 80
from k4c_driver import chamber_S
s3 = 1/mp.sqrt(3)
NC = [[s3,s3,s3],[s3,s3,-s3],[s3,-s3,s3],[-s3,s3,s3]]
J = json.load(open("research/08_theory/K4b_cubo_finish/tetra_Sdk_final.json"))
cert = [mpf(s) for s in J["Sdk"]]
for HI2 in [mpf(1)/64, mpf(1)/96]:
    t0=time.time()
    SC, st = chamber_S(NC, 4, mpf(1)/32, mpf(1)/32, mpf('4.5'), mpf('4.5'), HI2=HI2)
    el = time.time()-t0
    print(f"inner HI2=1/{int(1/HI2)} ({el:.0f}s) stats={st}")
    for k in range(5):
        print(f"   [d^{k}] err = {mp.nstr(abs(SC[k]-cert[k]),3)}", flush=True)
