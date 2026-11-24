
import os
ROOT = '<project>'
T = os.path.join(ROOT, 'research/08_theory')
print("===== deg4_residue_polynomials.json =====")
print(open(os.path.join(T, 'TOWN1_tower_large_n/claude_checks/deg4_residue_polynomials.json')).read())
print("\n===== WT5 kernels_deg4.json =====")
print(open(os.path.join(T, 'WT5_weight5_all_k/kernels_deg4.json')).read()[:3000])
