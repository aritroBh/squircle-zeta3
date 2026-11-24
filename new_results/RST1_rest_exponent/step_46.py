import os, subprocess, time
os.chdir('<project>/research/08_theory/RST1_rest_exponent')
t0=time.monotonic()
r = subprocess.run(['bash','-c','cd <project>/research/08_theory/RST1_rest_exponent && nice <project>/.venv/bin/python verify_rst1.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout)
print(r.stderr[-1000:] if r.stderr else '')
ok = "ALL RST1 CHECKS PASSED" in r.stdout
open('progress.log','a').write(f"verify_rst1.py full rerun: {'PASS' if ok else 'FAIL'} ({time.monotonic()-t0:.0f}s)\n")
