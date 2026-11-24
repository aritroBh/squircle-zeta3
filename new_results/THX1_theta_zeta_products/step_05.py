
import os
p = '<project>/research/08_theory/K3_fourier_all_orders'
for root, dirs, files in os.walk(p):
    for f in sorted(files):
        fp = os.path.join(root, f)
        print(fp, os.path.getsize(fp))
