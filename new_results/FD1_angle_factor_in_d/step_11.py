t2='<project>/research/08_theory/TOW2_tower_all_depths'
import re
txt=open(t2+'/proof_all_depths.md').read()
print(len(txt))
for m in re.finditer(r'^#+ .*', txt, re.M):
    print(m.start(), m.group(0))
