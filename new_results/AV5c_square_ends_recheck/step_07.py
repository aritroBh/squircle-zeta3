import os, glob
# find the paper
for p in glob.glob('../../../paper/ieee/*.tex'):
    print(p, os.path.getsize(p))
