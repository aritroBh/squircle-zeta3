import os, sys, time
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver, trunc_engine as te

D, K = 8, 38
outpath = f'data/towers_D{D}_K{K}.pkl'
t0 = time.time()
H = driver.harvest(D, K, outpath)
print(f"harvest done: {time.time()-t0:.1f}s, keys={len(H)}")
with open('progress.log','a') as f:
    f.write(f"step_06 harvest D={D} K={K}: {time.time()-t0:.0f}s, {len(H)} keys -> {outpath}\n")
