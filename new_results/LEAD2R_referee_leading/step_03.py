import os
for d in ['../LEAD2_leading_proof', '../LEAD1_leading_genfun', '../WTS_weight_k_structure', '../TOWN1_tower_large_n', '../TWOEND1_weight_poles']:
    print("="*60); print(d)
    for root, dirs, files in os.walk(d):
        for f in sorted(files):
            p = os.path.join(root, f)
            print("   ", p, os.path.getsize(p))
