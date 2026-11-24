import json
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
d9=json.load(open(fd1+'/Dtop_d9.json'))
print(json.dumps(d9, indent=1)[:2000])
d=json.load(open(fd1+'/Dtop.json'))
print(json.dumps({k:d[k] for k in list(d)[:3]}, indent=1))
