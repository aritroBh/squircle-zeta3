cwd = "<project>/research/08_theory/G2R_allN_proof_referee"
import os
for f in sorted(os.listdir(cwd)):
    p = os.path.join(cwd, f)
    print("=== ", f, os.path.getsize(p))
    print(open(p).read()[:2000])
