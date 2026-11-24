
base = '<project>/research'
import os
p = f'{base}/10_frontier/E1_edge_law'
for root, dirs, files in os.walk(p):
    dirs[:] = [d for d in dirs if d != '__pycache__']
    depth = root.replace(p,'').count(os.sep)
    print('  '*depth + os.path.basename(root) + '/')
    for f in sorted(files):
        print('  '*(depth+1) + f)
