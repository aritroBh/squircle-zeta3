import os, subprocess, sys, time
cwd = "<project>/research/08_theory/G2_M4alpha"
# status of flat_exact N520
out = open(cwd + "/flat_exact_N520_K4.out").read()
print(out[-400:])
print("N520 json exists:", os.path.exists(cwd + "/flat_exact_N520_K4.json"))
