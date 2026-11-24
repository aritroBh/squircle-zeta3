import os
base='<project>/research/08_theory'
p=os.path.join(base,'M4f_w5_extras')
print(sorted(os.listdir(p)))
r=os.path.join(p,'result.md')
if os.path.exists(r):
    print(open(r).read()[:5000])
