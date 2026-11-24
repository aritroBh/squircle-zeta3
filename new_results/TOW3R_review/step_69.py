import subprocess, os
base = '<project>/research/08_theory/'
d3 = base+'TOW3_tower_depth6/'
# find the actual Task C verification step (step_61 is just exploration). grep for 'cubic' or '40/40' or kern in steps 56-70
out = subprocess.run(['grep','-ln','kern_5\|kern_6\|ipos\|ineg'] + [d3+f'step_{n:02d}.py' for n in range(56,86)], capture_output=True, text=True)
print(out.stdout)
