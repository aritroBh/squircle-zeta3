base = '<project>/research/08_theory/'
import os
print(os.listdir(base+'WTS_weight_k_structure/theory/'))
t = open(base+'WTS_weight_k_structure/theory/result.md').read()
print(len(t))
import re
for m in re.finditer(r'^#+ .*', t, re.M): print(m.group(0))
