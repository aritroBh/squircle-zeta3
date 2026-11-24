import os
base = '<project>'
w6 = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'w6')
with open(os.path.join(w6, 'identify_w6.md')) as f:
    print(f.read()[:7000])
