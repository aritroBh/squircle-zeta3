import os
base = '<project>/research/08_theory'
src = open(os.path.join(base,'HEXC2R_referee_certificate','reflib.py')).read()
print(src[14000:19500])
