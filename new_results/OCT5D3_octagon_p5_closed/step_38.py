import subprocess, sys, os, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
r = subprocess.run([sys.executable, "gengine2.py", "1450", "560"], capture_output=True, text=True)
print(r.stdout.strip(), r.stderr[-500:] if r.returncode else "")
