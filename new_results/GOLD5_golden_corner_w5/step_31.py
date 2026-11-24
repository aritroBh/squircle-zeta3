import subprocess
files = ['../G2b_point_depth_M4/step_%02d.py' % i for i in range(1, 30)]
out = subprocess.run(['grep', '-ln', 'class\|P5\|Li5-span\|li5'] + files, capture_output=True, text=True)
print(out.stdout)
