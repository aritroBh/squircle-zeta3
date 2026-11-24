import os, json
base='<project>'
m4c=base+'/research/08_theory/M4c_w5_closed_form'
print([f for f in sorted(os.listdir(m4c)) if 'W5' in f or 'w5' in f][:40])
