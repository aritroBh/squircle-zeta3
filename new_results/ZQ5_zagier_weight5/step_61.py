import subprocess, sys, os
r = subprocess.run(['cat','progress.log'],capture_output=True,text=True)
print(r.stdout[-3000:])
print("=== out/d12:", sorted(os.listdir('out/d12')))
print("=== out/d5:", sorted(os.listdir('out/d5')))
print("=== out/d8:", sorted(os.listdir('out/d8')))
