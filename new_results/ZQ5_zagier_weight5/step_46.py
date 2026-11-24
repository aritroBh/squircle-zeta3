import sys, subprocess
src = r'''
import os
p = "../A3_W4_proof/w4_symbol.py"
print(os.path.exists(p))
t = open(p).read()
print("length:", len(t))
print(t[:9000])
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',src],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-500:])
