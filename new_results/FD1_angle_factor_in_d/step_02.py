import os
base='<project>'
a2=os.path.join(base,'research/08_theory/ALPHA2_angle_leading_law')
print(os.path.exists(a2))
for f in sorted(os.listdir(a2)):
    p=os.path.join(a2,f)
    print(f, os.path.getsize(p))
fd1=os.path.join(base,'research/08_theory/FD1_angle_factor_in_d')
print('FD1 exists:', os.path.exists(fd1))
if os.path.exists(fd1):
    print(os.listdir(fd1))
