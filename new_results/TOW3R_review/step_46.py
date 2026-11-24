base = '<project>/research/08_theory/'
t = open(base+'TOW1_tower_law/proof_tower.md').read()
print(len(t))
import re
for m in re.finditer(r'^#+ .*', t, re.M): print(m.group(0))
