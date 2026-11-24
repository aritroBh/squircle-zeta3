import os, shutil, json
l1='<project>/research/08_theory/LEAD1_leading_genfun'
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
print([f for f in os.listdir(l1) if 'top' in f.lower() or 'Dtop' in f])
for f in ['Dtop.json','Dtop_d9.json']:
    if os.path.exists(l1+'/'+f):
        shutil.copy(l1+'/'+f, fd1+'/'+f)
        d=json.load(open(fd1+'/'+f))
        print(f, type(d), list(d.keys())[:10] if isinstance(d,dict) else d[:5])
