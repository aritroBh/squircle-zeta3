import os, sys, time, shutil
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver

shutil.copy('data/towers_D8_K48.pkl', 'data/towers_D8_K52.pkl')
t0 = time.time()
H = driver.harvest(8, 52, 'data/towers_D8_K52.pkl')
print(f"extend to K=52: {time.time()-t0:.1f}s, keys={len(H)}")
with open('progress.log','a') as f:
    f.write(f"step_17 extend harvest D=8 to K=52 (+{time.time()-t0:.0f}s)\n")
