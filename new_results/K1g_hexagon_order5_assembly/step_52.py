import os
base = '<project>/research/08_theory'
e5c = os.path.join(base, 'E5c_nonlocal_z5')
with open(os.path.join(e5c,'symeng.py')) as f:
    c = f.read()
# Print the edge_params and log-related parts
import re
for i,ln in enumerate(c.split('\n')):
    if any(s in ln for s in ['def edge_params','log','nu','m_j','unit','conv','A0','LA','norm','|m','mag']):
        print(f"{i:3d}: {ln.rstrip()[:150]}")
