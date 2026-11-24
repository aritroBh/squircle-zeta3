import os, json, glob
target = '<project>/research/08_theory/COR8B_depth_two_two_terms'
# find hpl_fit file
for root, dirs, files in os.walk('<project>'):
    for f in files:
        if 'hpl_fit_k3' in f:
            print(os.path.join(root, f))
    # limit depth
    if root.count(os.sep) > 6:
        dirs[:] = []
