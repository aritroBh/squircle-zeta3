import os
cwd = "<project>/research/08_theory/G2_M4alpha"
print(sorted(os.listdir(cwd)))
print("---- current step_05.py head ----")
print(open(cwd + "/step_05.py").read()[:200])
