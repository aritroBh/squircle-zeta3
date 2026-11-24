import os
base = '<project>/research/08_theory'
e5c = os.path.join(base, 'E5c_nonlocal_z5')
with open(os.path.join(e5c,'symeng.py')) as f:
    c = f.read()
lines = c.split('\n')
for i in range(14, 60):
    print(f"{i:3d}: {lines[i].rstrip()}")
