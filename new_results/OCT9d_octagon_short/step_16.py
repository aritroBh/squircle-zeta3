import subprocess, time, sys
base = "<project>/research/08_theory/"
t0 = time.time()
r = subprocess.run([sys.executable, base+"OCT9c_octagon_close/verify_oct9c.py"], capture_output=True, text=True, cwd=base+"OCT9c_octagon_close")
print("elapsed", time.time()-t0)
print(r.stdout)
print(r.stderr[-2000:])
