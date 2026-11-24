import os, glob
base = os.getcwd()
# go up to find research dirs
for p in sorted(glob.glob('../*')):
    print(" ", p)
