import os
base='<project>/research/08_theory/M4c_w5_closed_form'
r=open(os.path.join(base,'result.md')).read()
print(r[:9000])
