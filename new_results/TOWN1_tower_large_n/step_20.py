import os, sys, time
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
TOW1 = os.path.join(base, 'research/08_theory/TOW1_tower_law')
sys.path.insert(0, TOW1); sys.path.insert(0, WTS)
import trunc_engine as te
import towlib

# Validate at D=8, k<=14 (both rho and tower)
S2 = te.stirling_table(16)
D = 8
t0 = time.time()
bad = 0; tot = 0
for k in range(2, 15):
    rho_k = te.rho_level(k, D, S2)
    tower_k = te.build_tower_level(k, D, rho_k)
    for n in range(2, k+1):
        for key in te.keys_of_degree(n):
            p, q = key
            for delta in range(0, min(D, 2*n-4)+1):
                ref = towlib.u(p, q, k, delta)
                if ref != tower_k[key][delta]:
                    bad += 1
                tot += 1
print(f"D=8 tower check: {tot} coeffs k<=14, mismatches={bad}, time={time.time()-t0:.2f}s")
with open('progress.log','a') as f:
    f.write(f"step_05 D=8 tower validated vs towlib.u: {tot} coeffs k<=14, {bad} mismatch\n")

# timing for larger k at D=8
for k in [20, 26, 32]:
    t0 = time.time()
    rho_k = te.rho_level(k, D, te.stirling_table(k+1))
    t1 = time.time()
    tower_k = te.build_tower_level(k, D, rho_k)
    t2 = time.time()
    print(f"k={k}: rho={t1-t0:.2f}s tower={t2-t1:.2f}s  keys_deg_k={len(te.keys_of_degree(k))}")
