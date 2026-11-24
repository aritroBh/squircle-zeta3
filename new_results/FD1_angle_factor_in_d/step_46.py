a1='<project>/research/08_theory/ALPHA1_angle_family_across_k'
import os
print([f for f in os.listdir(a1) if f.startswith('step') or f.endswith('.py')])
print("=== step_07 ===")
print(open(a1+'/step_07.py').read()[:2500])
