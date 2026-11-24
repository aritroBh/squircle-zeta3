base = "<project>/research/08_theory"
import os
for f in sorted(os.listdir(f"{base}/OCT5_octagon_exact")):
    print(f)
print("----- brief:")
with open(f"{base}/OCT5_octagon_exact.brief") as f:
    print(f.read()[:4000])
