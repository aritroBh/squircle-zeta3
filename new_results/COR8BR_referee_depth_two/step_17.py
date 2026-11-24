import os
target = '<project>/research/08_theory/COR8B_depth_two_two_terms'
with open(os.path.join(target, 'progress.log')) as f:
    print(f.read())
