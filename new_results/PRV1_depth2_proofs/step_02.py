
import os
base='research/08_theory/DEC1_decagon'
for root,dirs,files in os.walk(base):
    print(root, files)
