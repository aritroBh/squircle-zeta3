import os, subprocess
base = '<project>/research/08_theory/'
print(os.listdir(base))
print(os.listdir(base + 'TOW3_tower_depth6/'))
work = base + 'TOW3R_review/'
os.makedirs(work, exist_ok=True)
print(os.listdir(work) if os.path.exists(work) else 'no dir')
