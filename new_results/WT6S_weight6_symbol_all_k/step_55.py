
import subprocess, os
ROOT = '<project>'
# grep for the c(6) constant digits across research
r = subprocess.run(['grep', '-rl', '335687558224020357', os.path.join(ROOT, 'research')],
                   capture_output=True, text=True, timeout=120)
print(r.stdout)
print("stderr", r.stderr[:500])
