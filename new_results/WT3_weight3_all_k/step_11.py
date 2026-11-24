
import os
BASE = '<project>'
WT2 = os.path.join(BASE, 'research/08_theory/WT2_weight2_all_k')
print("===== genintegrand.py =====")
print(open(os.path.join(WT2,'genintegrand.py')).read())
print("\n===== buildintegrand.py =====")
print(open(os.path.join(WT2,'buildintegrand.py')).read())
