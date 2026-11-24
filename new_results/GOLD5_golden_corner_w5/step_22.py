import subprocess
out = subprocess.run(['grep', '-l', 'sols_golden\|golden_sols\|unit_eq\|x+y\|x + y', '../COR9_golden_corner_closed_forms/'+f for f in ['step_%02d.py'%i for i in range(1,30)]], capture_output=True, text=True)
print(out.stdout, out.stderr)
out = subprocess.run(['grep', '-rn', 'phi', '../COR9_golden_corner_closed_forms/step_01.py'], capture_output=True, text=True)
print(out.stdout[:2000])
