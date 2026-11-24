base = '<project>/research/08_theory'
tn = f'{base}/TOWN1_tower_large_n'
print(open(f'{tn}/driver.py').read()[:4000])
import os
print(os.listdir(f'{tn}/data'))
