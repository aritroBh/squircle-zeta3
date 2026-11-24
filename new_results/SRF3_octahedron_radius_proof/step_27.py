
import subprocess, time, os
PY = "<project>/.venv/bin/python"
h1d = "<project>/research/10_frontier/H1_octahedron_end"
for K in [20, 40]:
    t0=time.time()
    r = subprocess.run(["nice","-n","10",PY,"05_high_order_numeric.py","3",str(K)],capture_output=True,text=True,cwd=h1d,timeout=2400)
    print("K=",K," time: %.1fs"%(time.time()-t0), r.stdout.splitlines()[-1] if r.stdout else r.stderr[-300:])
