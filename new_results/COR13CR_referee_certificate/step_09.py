import os
base = '<project>/research/08_theory'
src = open(os.path.join(base,'HEXC2R_referee_certificate','reflib.py')).read()
# find distribution definition
i = src.find('def distribution')
print(src[i:i+4200])
