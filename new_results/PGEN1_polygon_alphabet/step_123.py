# step_94: run verify_pgen1.py
import subprocess, time
t0 = time.time()
r = subprocess.run(["nice", "-n", "10", "<project>/.venv/bin/python", "verify_pgen1.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/PGEN1_polygon_alphabet")
print("returncode:", r.returncode, " wall:", round(time.time()-t0,1), "s")
print(r.stdout[-6000:])
if r.returncode != 0:
    print("STDERR:", r.stderr[-2000:])
