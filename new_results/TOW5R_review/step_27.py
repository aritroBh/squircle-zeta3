import os
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
town1 = '<project>/research/08_theory/TOWN1_tower_large_n'
tw = sorted(f for f in os.listdir(wts) if f.startswith('tower_') or f.startswith('binteg_') or f.startswith('kern_'))
print('WTS tower_k:', sorted(int(f.split('_')[1].split('.')[0]) for f in tw if f.startswith('tower_')))
print('WTS binteg_d:', sorted(int(f.split('_')[1].split('.')[0]) for f in tw if f.startswith('binteg_')))
print('WTS kern_k:', sorted(int(f.split('_')[1].split('.')[0]) for f in tw if f.startswith('kern_')))
print()
print('TOWN1 files:', sorted(os.listdir(town1))[:40])
if os.path.exists(town1+'/data'):
    print('TOWN1 data:', sorted(os.listdir(town1+'/data')))
