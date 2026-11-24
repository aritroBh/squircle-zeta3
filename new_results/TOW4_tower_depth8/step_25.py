import os, pickle, time
import sympy as sp
base = '<project>/research/08_theory'
town1 = base + '/TOWN1_tower_large_n'
print(sorted(os.listdir(town1)))
with open(town1 + '/driver.py') as f:
    drv = f.read()
# show format documentation
i0 = drv.find('pkl')
print(drv[:3000])
