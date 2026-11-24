
import os
BASE = '<project>'
res = open(os.path.join(BASE, 'research/08_theory/WT2b_stabilization/result.md')).read()
print(res[:10000])
