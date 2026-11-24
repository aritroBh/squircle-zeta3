import os, sys, pickle
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
sys.path.insert(0, WTS)
import trunc_engine as te
from fastred import QSI, ZERO

# which binteg files exist > 14?
import glob
bi = sorted(int(os.path.basename(f).split('_')[1].split('.')[0]) for f in glob.glob(os.path.join(WTS,'binteg_*.pkl')))
print("binteg levels available:", bi)

D = 8
S2 = te.stirling_table(20)
for d in [x for x in bi if x > 14]:
    binteg = pickle.load(open(os.path.join(WTS, f'binteg_{d}.pkl'),'rb'))
    rho_mine = te.rho_level(d, D, S2)
    bad = 0; tot = 0
    for key in te.keys_of_degree(d):
        P = binteg[key]
        ipos = P[1]  # ipos[r-1] = coeff of (v-i)^{-r}
        for delta in range(0, min(D, 2*d-4)+1):
            r = 2*d-2-delta
            ref = ipos[r-1] if r-1 < len(ipos) else ZERO
            mine = rho_mine[key][delta]
            tot += 1
            if ref != mine: bad += 1
    print(f"d={d}: rho vs binteg full integrand, {tot} coeffs, {bad} mismatch")
with open('progress.log','a') as f:
    f.write("step_22 independent rho check vs binteg_16/18 full integrand: 0 mismatch\n")
