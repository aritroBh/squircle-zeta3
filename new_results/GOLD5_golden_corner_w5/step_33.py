import subprocess
files = ['../G2b_point_depth_M4/step_%d.py' % i for i in range(1, 117)]
out = subprocess.run(['grep', '-ln', 'P5\|li5_basis\|Li5 basis\|class solve\|classify'] + files, capture_output=True, text=True)
print(out.stdout)
