import os, glob
base='<project>/research/08_theory'
for d in glob.glob(base+'/TOW2*')+glob.glob(base+'/LEAD2*'):
    print(d, '->', os.listdir(d)[:40] if os.path.isdir(d) else 'file')
