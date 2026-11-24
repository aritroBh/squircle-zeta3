import os
base='research/10_frontier/E1_edge_law'
src=open(os.path.join(base,'s04_exact_local.py')).read()
print(src[1500:2900])
