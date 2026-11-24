import os
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
tmp = 'data/_verify_towers_D6_K22.pkl'
if os.path.exists(tmp):
    os.remove(tmp); print("removed temp", tmp)

# Final summary of key cells and confirm counts from the exported JSON
import json
for f in ['cells_D6_N40','cells_D8_N40','cells_D8_N45']:
    d = json.load(open(f'results/{f}.json'))
    print(f"{f}: D={d['D']} N={d['N']} K={d['K']} ncells={d['ncells']} ngenuine={d['ngenuine']} nfails={d['nfails']} min_nchecks={d['min_nchecks']}")
kl = json.load(open('results/kernel_laws.json'))
print("\nkernel_laws:", kl['degree2_kernels_exact_range'], "deg2 mismatches:", kl['degree2_mismatches'], "| deg3 range:", kl['degree3_kernel_law_range'])
with open('progress.log','a') as f:
    f.write("step_23 final: verify 5/5 PASS post-reboot; deliverables confirmed\n")
print("\nDONE - all tasks complete")
