import os, glob, re
base='<project>/research/08_theory'
cands = glob.glob(base+'/ALPHA1*')+glob.glob(base+'/LEAD2*')+glob.glob(base+'/LEAD1*')
for d in cands:
    print(d)
