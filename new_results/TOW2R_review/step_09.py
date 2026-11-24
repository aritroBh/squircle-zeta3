import os
base = '<project>/research/08_theory'
tow2 = os.path.join(base, 'TOW2_tower_all_depths')
workdir = os.path.join(base, 'TOW2R_review')
os.chdir(workdir)
with open('progress.log', 'a') as f:
    f.write('step_01: read TOW2 result.md + proof_all_depths.md fully; read tow2lib.py, prove_S12/S34/tops/adj34/E3b34.py\n')

# Read TOW1 proof (base). Check size first
tow1 = os.path.join(base, 'TOW1_tower_law')
print(os.listdir(tow1))
