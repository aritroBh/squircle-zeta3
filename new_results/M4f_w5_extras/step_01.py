import os, subprocess
# Find our bearings
for root, dirs, files in os.walk('<project>'):
    depth = root.replace('<project>','').count(os.sep)
    if depth <= 3 and ('M4' in root or root.endswith('Research')):
        print(root, len(files), 'files')
