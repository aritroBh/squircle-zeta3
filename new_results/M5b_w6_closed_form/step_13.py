import os, json
base='<project>/research/08_theory'
r=open(os.path.join(base,'M4f_w5_extras/result.md')).read()
print(r[:12000])
