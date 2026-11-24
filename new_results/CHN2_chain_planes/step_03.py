
import os
print("planner dir:")
for f in sorted(os.listdir('planner')):
    print("  ", f)
print()
# Find CHN1
base = '<project>/research/08_theory'
for d in sorted(os.listdir(base)):
    print(d)
