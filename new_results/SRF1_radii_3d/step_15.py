
import os
f = '<project>/research/10_frontier'
for root, dirs, files in os.walk(f):
    print(root)
    for x in files: print('   ', x)
