import os
base = '<project>/research/08_theory'
src = open(os.path.join(base,'HEXC2R_referee_certificate','steps','step_06.py')).read()
i = src.find('# (d)')
print(src[i:i+2200] if i>=0 else src[-2200:])
