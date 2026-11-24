import os
base = '<project>'
w4 = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'w4')
with open(os.path.join(w4, 'identify_w4.md')) as f:
    print(f.read()[:6000])
