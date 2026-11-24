import os
base = "<project>/research/"
cwd = "<project>/research/08_theory/G2R_allN_proof_referee"
print(os.getcwd())
for root, dirs, files in os.walk(os.path.join(base, "08_theory/G2_M4alpha")):
    for f in sorted(files):
        p = os.path.join(root, f)
        print(p, os.path.getsize(p))
