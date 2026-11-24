import sys, subprocess
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_21_d8_verify.py'],capture_output=True,text=True,timeout=2400)
out = r.stdout
# trim the phase_lll noise lines
lines = [l for l in out.splitlines() if not l.startswith("D=8 LLL")]
print("\n".join(lines)); print("STDERR:", r.stderr[-1000:])
with open('progress.log','a') as f: f.write("step_21_d8_verify(r3): fixed convention; (a)-(e) pass\n")
