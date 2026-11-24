import subprocess, os
os.chdir('<project>/research/08_theory/PGEN1R_referee_alphabet')
# Initialize progress log and look at PGEN1 steps 45-53 headers to find N3 definition
for i in range(45, 54):
    fn = f'../PGEN1_polygon_alphabet/step_{i:02d}.py'
    if not os.path.exists(fn): continue
    with open(fn) as f:
        head = f.read(600)
    print('='*30, fn)
    print(head)
