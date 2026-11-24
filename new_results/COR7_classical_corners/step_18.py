import os
base = '<project>/research/'
# find hpl_series.py anywhere
for root, dirs, files in os.walk(base):
    for f in files:
        if 'hpl_series' in f or 'hpl_fit' in f:
            print(os.path.join(root, f))
