t2='<project>/research/08_theory/TOW2_tower_all_depths'
txt=open(t2+'/result.md').read()
print(len(txt))
# find the section about top recursion / e(2s,s)
import re
idx=txt.find('e(2s')
print(txt[:200])
# print table of contents-ish
for m in re.finditer(r'^#+ .*', txt, re.M):
    print(m.group(0))
