base = '<project>/research/08_theory'
import os
for f in ['proof_all_depths.md', 'tow2lib.py']:
    p = f'{base}/TOW2_tower_all_depths/{f}'
    print(f, os.path.getsize(p), 'bytes')
# Also check TOW3 existing step_01
with open(f'{base}/TOW3_tower_depth6/step_01.py') as f:
    print(f.read())
