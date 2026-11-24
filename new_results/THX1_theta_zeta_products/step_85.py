
import subprocess, time
t=time.time()
r = subprocess.run(['<project>/.venv/bin/python', 'verify_thx1.py'],
                   capture_output=True, text=True, cwd='.')
print(r.stdout)
print("STDERR:", r.stderr[-2000:] if r.stderr else "")
print("wall time: %.1fs" % (time.time()-t))
