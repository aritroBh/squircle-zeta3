import os
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
print(sorted(os.listdir(theory)))
# Read genfun_d1_d3.md
with open(os.path.join(theory, 'genfun_d1_d3.md')) as f:
    print(f.read())
