import os
base = '../PGEN1R_referee_alphabet'
for p in sorted(os.listdir(base)):
    full = os.path.join(base, p)
    print(f"{p:40s} {os.path.getsize(full):>10d}")
print()
print(open('step_01.py').read())
