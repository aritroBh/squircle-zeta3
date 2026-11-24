
# Peek at K2 result.md section 3 to get the exact representation
txt = open('research/08_theory/K2_corner_decay/result.md').read()
print(len(txt))
# find section 3
import re
m = re.search(r'##\s*3.*?(?=\n##\s*4|\Z)', txt, re.S)
print(m.group(0)[:6000] if m else txt[:6000])
