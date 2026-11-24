import os
base = '<project>'
theory = os.path.join(base, 'research', '08_theory', 'WTS_weight_k_structure', 'theory')
with open(os.path.join(theory, 'run_kernel_k.py')) as f:
    src = f.read()
print(src[:4000])
