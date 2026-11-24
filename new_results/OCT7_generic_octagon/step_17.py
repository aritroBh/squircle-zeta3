import os
base="research/08_theory"
for d in ["E5_edge_law","K1_nonlocal_law","D15_polygon_p1"]:
    p=os.path.join(base,d)
    print("="*25, d, "="*25)
    print(sorted(os.listdir(p))[:40])
