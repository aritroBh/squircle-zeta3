import os
d = 'research/08_theory/DEC1_decagon'
for root, dirs, files in os.walk(d):
    print(root)
    for f in files: print('   ', f)
