import os
base='research/10_frontier/E1_edge_law'
print(os.listdir(base))
print(open(os.path.join(base,'result.md')).read()[:8000])
