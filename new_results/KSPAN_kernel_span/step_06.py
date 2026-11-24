import os, sys

wt2_path = "../WT2_weight2_all_k"
if os.path.exists(wt2_path):
    print(f"WT2_weight2_all_k contents:")
    for item in sorted(os.listdir(wt2_path)):
        print(f"  {item}")
else:
    print(f"WT2_weight2_all_k doesn't exist at {wt2_path}")