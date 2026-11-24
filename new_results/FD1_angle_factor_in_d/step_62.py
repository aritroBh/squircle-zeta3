t2='<project>/research/08_theory/TOW2_tower_all_depths'
src=open(t2+'/tow2lib.py').read()
print(len(src))
import re
for m in re.finditer(r'def \w+', src): print(m.group(0))
# print sections about zigzag
i=src.find('zigzag')
print(src[max(0,i-200):i+3000])
