import os
base = '<project>/research/08_theory'
src = open(os.path.join(base,'COR13C_third_corner_certificate','cor13c_common.py')).read()
i = src.find('def cr(')
print(src[i:i+3200])
