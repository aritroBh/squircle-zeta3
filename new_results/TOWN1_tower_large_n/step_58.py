import subprocess, os, pickle, time
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
PY = base + '/.venv/bin/python'

# data integrity check
for f in ['data/towers_D8_K38.pkl','data/towers_D8_K48.pkl','data/towers_D8_K52.pkl']:
    H = pickle.load(open(f,'rb'))
    nk = len([k for k in H if k != '_meta'])
    print(f"{f}: loads OK, {nk} keys, meta={H.get('_meta')}")

t0=time.time()
r = subprocess.run(['nice','-n','10',PY,'verify_town1.py'], capture_output=True, text=True, cwd=workdir)
print("\n--- verify_town1.py tail ---")
print('\n'.join(r.stdout.splitlines()[-8:]))
print("rc:", r.returncode, f"({time.time()-t0:.1f}s)")
