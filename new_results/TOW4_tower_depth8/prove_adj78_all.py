import subprocess, sys, time
PY = sys.executable
work = '<project>/research/08_theory/TOW4_tower_depth8'
cases = [(7,0),(7,1),(7,2),(7,3),(8,0),(8,1),(8,2),(8,3),(8,4)]
results = {}
for delta,s in cases:
    t0=time.time()
    r = subprocess.run([PY, work+'/prove_adj78_fast.py',str(delta),str(s)], capture_output=True, text=True)
    ok = 'True' in r.stdout
    results[(delta,s)] = ok
    line = [l for l in r.stdout.splitlines() if 'adjacent' in l]
    print(f"({delta},{s}): {'PASS' if ok else 'FAIL'} {line[-1] if line else r.stderr[-200:]} wall {round(time.time()-t0,1)}s", flush=True)
print("ALL:", all(results.values()), results, flush=True)
