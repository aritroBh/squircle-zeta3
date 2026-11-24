import subprocess, sys, os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
r = subprocess.run([sys.executable, "oct5d3_basis.py"], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else "OK")
