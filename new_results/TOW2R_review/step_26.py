import os
base = '<project>/research/08_theory'
wts = os.path.join(base, 'WTS_weight_k_structure', 'theory')
print(sorted(os.listdir(wts))[:40])
town1 = os.path.join(base, 'TOWN1_tower_large_n')
print(os.path.exists(town1), os.listdir(town1) if os.path.exists(town1) else '')
if os.path.exists(os.path.join(town1, 'data')):
    print(sorted(os.listdir(os.path.join(town1,'data'))))
