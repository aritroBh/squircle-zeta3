code = open('steps/step_07.py').read()
old = """        if sums is None: sums = [mp.mpf(0)]*6
        for j in js:"""
new = """        if sums is None: sums = [mp.mpf(0)]*6
        else: sums = [s/2 for s in sums]   # node reuse: old nodes carry half weight at new du
        for j in js:"""
assert old in code
code = code.replace(old, new)
open('steps/step_07.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_07.py','80'],capture_output=True,text=True,timeout=5400)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-2000:])
