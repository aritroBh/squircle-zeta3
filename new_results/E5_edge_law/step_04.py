import os
base='research/10_frontier/E1_edge_law'
src=open(os.path.join(base,'s04_exact_local.py')).read()
# print the middle part that was truncated
print(src[2800:5200])
