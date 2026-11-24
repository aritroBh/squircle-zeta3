cwd = "<project>/research/08_theory/G2_M4alpha"
import os, time
print(open(cwd + "/flat_exact_N400_K4.out").read()[-700:])
print("json size", os.path.getsize(cwd + "/flat_exact_N400_K4.json"))
