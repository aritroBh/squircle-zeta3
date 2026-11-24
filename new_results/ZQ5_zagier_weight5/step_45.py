import sys, subprocess
src = r'''
import os
d = "../M4i_w5_noL"
for fn in ["w5_symbol.py", "lrel_verify.py", "nice_symbol.json"]:
    p = os.path.join(d, fn)
    print("="*30, fn, os.path.exists(p))
    t = open(p).read()
    print(t[:5000])
'''
r = subprocess.run(['nice','-n','10',sys.executable,'-c',src],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-500:])
