import os, sys, subprocess
src = r'''
import os
print("cwd:", os.getcwd())
print("== parent (08_theory) ==")
for f in sorted(os.listdir("..")):
    p = os.path.join("..", f)
    print(("D " if os.path.isdir(p) else "F "), f)
print("== grandparent (research) ==")
for f in sorted(os.listdir("../.."))[:60]:
    p = os.path.join("../..", f)
    print(("D " if os.path.isdir(p) else "F "), f)
for cand in ["../M4i_w5_noL", "../M4i_w5_noL/L_relation.json"]:
    print(cand, os.path.exists(cand))
if os.path.exists("step_01.py"):
    print("== existing step_01.py (first 2000 chars) ==")
    print(open("step_01.py").read()[:2000])
'''
open('steps/step_02_ls.py','w').write(src)
r = subprocess.run(['nice','-n','10',sys.executable,'steps/step_02_ls.py'],capture_output=True,text=True,timeout=300)
print(r.stdout); print("STDERR:", r.stderr[-800:])
with open('progress.log','a') as f: f.write("step_02_ls: parent listing done\n")
