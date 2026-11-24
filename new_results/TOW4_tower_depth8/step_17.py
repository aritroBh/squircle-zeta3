import os
base = '<project>/research/08_theory'
print("TOW4 dir:", sorted(os.listdir(base + '/TOW4_tower_depth8')))
print("TOW2 dir:", sorted(os.listdir(base + '/TOW2_tower_all_depths'))[:40])
