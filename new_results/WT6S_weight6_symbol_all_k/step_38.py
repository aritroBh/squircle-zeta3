
import os, json
ROOT = '<project>'
W5 = os.path.join(ROOT, 'research/08_theory/WT5_weight5_all_k')
print(open(os.path.join(W5, 'pairings.json')).read())
