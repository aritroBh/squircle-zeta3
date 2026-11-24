import subprocess, os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_43.py'],
                   capture_output=True, text=True, timeout=1750)
print(r.stdout); print("STDERR:", r.stderr[-500:] if r.returncode!=0 else "(ok)")
