
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
