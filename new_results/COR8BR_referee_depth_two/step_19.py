import os
target = '<project>/research/08_theory/COR8B_depth_two_two_terms/steps'
for f in sorted(os.listdir(target)):
    print('='*40, f)
    with open(os.path.join(target, f)) as fh:
        print(fh.read())
