import os, subprocess, sys
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
for f in os.listdir("."):
    if f.startswith("gtab2_"): os.remove(f)
r = subprocess.run([sys.executable, "gengine2p.py", "430", "1440", "12"], capture_output=True, text=True)
print(r.stdout.strip(), r.stderr[-800:] if r.returncode else "")
