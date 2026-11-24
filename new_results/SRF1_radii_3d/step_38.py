
import os
p7 = '<project>/research/07_push'
for root, dirs, files in os.walk(p7):
    print(root)
    for x in files: print('   ', x)
