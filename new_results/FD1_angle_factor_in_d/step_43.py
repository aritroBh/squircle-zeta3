import os, glob
base='<project>'
hits=[]
for root,dirs,files in os.walk(base+'/research'):
    for f in files:
        if 'kern' in f and f.endswith('.pkl'):
            hits.append(os.path.join(root,f))
for h in hits[:60]: print(h, os.path.getsize(h))
