import os, sys, time
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
sys.path.insert(0, WTS)
import trunc_engine as te

D = 8
S2 = te.stirling_table(48)
for k in [36, 40, 46]:
    t0 = time.time()
    rho_k = te.rho_level(k, D, S2)
    t1 = time.time()
    tower_k = te.build_tower_level(k, D, rho_k)   # full tower to degree 2
    t2 = time.time()
    print(f"k={k}: rho={t1-t0:.2f}s fulltower={t2-t1:.2f}s")
