code = open('steps/step_07.py').read()
code = code.replace("def tanh_sinh_edge(L, h, s0, tol, du0=mp.mpf('0.09'), umax=mp.mpf('3.6'), levels=5):",
                    "def tanh_sinh_edge(L, h, s0, tol, du0=mp.mpf('0.09'), umax=mp.mpf('5.0'), levels=5):")
code = code.replace("""        if prev is not None:
            diff = max(abs(sums[k]-prev[k]) for k in range(6))
            if diff < tol: break
        prev = list(sums)""",
"""        if prev is not None:
            diff = max(abs(sums[k]-prev[k]) for k in range(6))
            print(f"    level {lev}: E0 err={mp.nstr(abs(sums[0]-L),3)} level-diff={mp.nstr(diff,3)}")
            if diff < tol: break
        prev = list(sums)""")
open('steps/step_07.py','w').write(code)
import subprocess, time
t0=time.time()
r=subprocess.run(['<project>/.venv/bin/python','steps/step_07.py','80'],capture_output=True,text=True,timeout=7200)
print("wall: %.1f s"%(time.time()-t0))
print(r.stdout); print("ERR:", r.stderr[-1500:])
