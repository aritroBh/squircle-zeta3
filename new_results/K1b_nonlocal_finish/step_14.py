import os
base = "<project>/research/10_frontier/E1_edge_law/"
for root, dirs, files in os.walk(base):
    depth = root[len(base):].count(os.sep)
    if depth < 2:
        print(root, sorted(files)[:25])
