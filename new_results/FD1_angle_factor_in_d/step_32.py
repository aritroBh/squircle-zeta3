import glob, os, re
base='<project>/research/08_theory'
t1=glob.glob(base+'/TOW1*')
print(t1)
for d in t1:
    if os.path.isdir(d):
        print(d, [f for f in os.listdir(d) if f.endswith('.md') or f.endswith('.pkl')][:20])
