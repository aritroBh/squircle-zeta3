import os
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
with open(os.path.join(theory, 'result.md')) as f:
    print(f.read())
