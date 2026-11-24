
import os
os.chdir("research/08_theory/D13f_sweep")
print("=== FILES IN D13f_sweep ===")
for fn in sorted(os.listdir(".")):
    print(fn, os.path.getsize(fn))
