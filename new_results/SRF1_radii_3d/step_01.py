
import os, subprocess
# Explore the research directory structure
base = '<project>'
for root, dirs, files in os.walk(base):
    depth = root.replace(base, '').count(os.sep)
    if depth <= 2:
        print(root)
        for f in files[:20]:
            print('   ', f)
    if depth > 2:
        dirs[:] = []
