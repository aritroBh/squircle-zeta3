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

S2 = te.stirling_table(16)
D = 6
t0 = time.time()
bad = []
tot = 0
for k in range(2, 15):
    rho_k = te.rho_level(k, D, S2)
    tower_k = te.build_tower_level(k, D, rho_k)
    # compare with recorded tower for all keys, depths <= D
    for key in te.keys_of_degree_all(k) if hasattr(te,'keys_of_degree_all') else []:
        pass
    for n in range(2, k+1):
        for key in te.keys_of_degree(n):
            p, q = key
            for delta in range(0, min(D, 2*n-4)+1):
                ref = towlib.u(p, q, k, delta)
                mine = tower_k[key][delta]
                tot += 1
                if ref != mine:
                    bad.append((k, key, delta, ref, mine))
print(f"tower check: {tot} coefficients (k=2..14, all keys, depths<=min({D},2n-4)), mismatches = {len(bad)}, time={time.time()-t0:.2f}s")
for b in bad[:10]:
    print("BAD", b)
with open('progress.log','a') as f:
    f.write(f"step_03 tower recursion validated vs towlib.u: {tot} coeffs k<=14 all keys depths<=D, {len(bad)} mismatch\n")
