
import os, subprocess
# find hpl_fit file
for root, dirs, files in os.walk('../..'):
    for f in files:
        if 'hpl_fit' in f.lower():
            print(os.path.join(root, f))
