base='<project>/research/08_theory'
W=base+'/CHN5_p5_w12_exact'
src=open(W+'/mzv_exact.py').read()
# print the BASIS section
import re
i=src.find('BASIS')
print(src[i-200:i+3000])
