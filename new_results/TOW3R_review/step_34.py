import subprocess, os
base = '<project>/research/08_theory/'
out = subprocess.run(['find', base, '-name', 'genintegrand.py', '-not', '-path', '*/node_modules/*'], capture_output=True, text=True)
print(out.stdout)
out2 = subprocess.run(['find', base, '-name', 'buildintegrand.py'], capture_output=True, text=True)
print(out2.stdout[:2000])
