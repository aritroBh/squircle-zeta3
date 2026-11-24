
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
src = open(G2B + "step_07.py").read()
# drop the broken triple search entirely
i = src.find("    # also record best triple")
j = src.find("    json.dump({")
src = src[:i] + src[j:]
open(G2B + "step_07.py", "w").write(src)
import subprocess, time
PY = "<project>/.venv/bin/python"
t0 = time.time()
r = subprocess.run(["nice", "-n", "19", PY, "step_07.py"], cwd=G2B, capture_output=True, text=True, timeout=1800)
print("elapsed", time.time() - t0, "rc", r.returncode)
print(r.stdout)
print("STDERR:", r.stderr[-1500:])
