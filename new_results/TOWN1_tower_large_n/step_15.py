import os, sys, time
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
TOW1 = os.path.join(base, 'research/08_theory/TOW1_tower_law')
sys.path.insert(0, TOW1); sys.path.insert(0, WTS)
import importlib
import trunc_engine as te
importlib.reload(te)
import towlib

S2 = te.stirling_table(16)
D = 6
t0 = time.time()
bad = []
tot = 0
for d in range(2, 15):
    rho_mine = te.rho_level(d, D, S2)
    for key in te.keys_of_degree(d):
        p, q = key
        for delta in range(0, min(D, 2*d-4)+1):
            ref = towlib.rho(key, delta)
            mine = rho_mine[key][delta]
            tot += 1
            if ref != mine:
                bad.append((key, delta, ref, mine))
print(f"rho check: {tot} coefficients (d=2..14, depths<=min({D},2d-4)), mismatches = {len(bad)}, time={time.time()-t0:.2f}s")
for b in bad[:10]:
    print("BAD", b)
with open('progress.log','a') as f:
    f.write(f"step_02 rho validated vs towlib.rho: {tot} coeffs d<=14 depths<=D, {len(bad)} mismatch\n")
