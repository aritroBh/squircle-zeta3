
import os
root='research'
hits=[]
for dp,dn,fn in os.walk(root):
    for f in fn:
        if 'D13' in f or 'd13' in f.lower():
            hits.append(os.path.join(dp,f))
for h in sorted(hits)[:40]: print(h)
