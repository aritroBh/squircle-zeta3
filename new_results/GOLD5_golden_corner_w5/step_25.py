import subprocess
files = ['../COR9_golden_corner_closed_forms/step_%02d.py' % i for i in range(1, 96)]
out = subprocess.run(['grep', '-ln', 'sunit\|s_unit\|Sunit\|solutions'] + files, capture_output=True, text=True)
print(out.stdout)
out = subprocess.run(['grep', '-n', 'gol\|golden', '../COR9_golden_corner_closed_forms/progress.log'], capture_output=True, text=True)
print(out.stdout[:3000])
