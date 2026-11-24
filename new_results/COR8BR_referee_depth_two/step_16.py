import os, json
target = '<project>/research/08_theory/COR8B_depth_two_two_terms'
with open(os.path.join(target, 'DT_A_residual_2kept.json')) as f:
    print(f.read())
print('='*50)
with open(os.path.join(target, 'classical_identity.txt')) as f:
    print(f.read()[:5000])
