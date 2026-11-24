import subprocess, os
base='<project>'
r=subprocess.run(['find', base+'/research', '-name', 'w4_coproduct.py'], capture_output=True, text=True)
print(r.stdout)
r=subprocess.run(['find', base, '-name', 'qsqrt2.py'], capture_output=True, text=True)
print(r.stdout)
