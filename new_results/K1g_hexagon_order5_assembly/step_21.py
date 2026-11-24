import os, json
base = '<project>/research/08_theory'
e5d = os.path.join(base, 'E5d_nonlocal_pslq')
print("=== E5d contents ===")
for d in sorted(os.listdir(e5d)):
    print(d)
