import os
base = '<project>'
p = os.path.join(base, 'research/08_theory/WT2_weight2_all_k')
with open(os.path.join(p, 'fastred.py')) as f:
    fastred = f.read()
print(fastred[9000:16000])
