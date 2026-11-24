import glob, json
G1 = "../../10_frontier/G1_zeta3_any_angle"
import os
print(os.path.exists(G1), os.getcwd())
for p in sorted(glob.glob(G1 + "/*")):
    print("  ", os.path.basename(p), os.path.getsize(p))
