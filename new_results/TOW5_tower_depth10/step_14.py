import os
town1 = '<project>/research/08_theory/TOWN1_tower_large_n'
for f in sorted(os.listdir(town1)):
    print(f, os.path.getsize(os.path.join(town1,f)))
print("="*80)
with open(os.path.join(town1,'driver.py')) as f:
    print(f.read())
