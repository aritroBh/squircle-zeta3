import os, sys, time, shutil, pickle
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver

# copy K38 harvest to K48 path and resume
shutil.copy('data/towers_D8_K38.pkl', 'data/towers_D8_K48.pkl')
t0 = time.time()
H = driver.harvest(8, 48, 'data/towers_D8_K48.pkl')
print(f"extend to K=48 done: {time.time()-t0:.1f}s, keys={len(H)}")
with open('progress.log','a') as f:
    f.write(f"step_07 extend harvest D=8 to K=48 (+{time.time()-t0:.0f}s)\n")
