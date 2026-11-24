import subprocess
files = ['../COR9_golden_corner_closed_forms/step_%02d.py' % i for i in range(1, 30)]
out = subprocess.run(['grep', '-l', 'unit'] + files, capture_output=True, text=True)
print(out.stdout)
out = subprocess.run(['grep', '-n', 'basis', '../COR9_golden_corner_closed_forms/step_04.py'], capture_output=True, text=True)
print(out.stdout[:3000])
