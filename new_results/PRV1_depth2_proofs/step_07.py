
import os
base='research/08_theory/DEC1_decagon'
txt = open(base+'/result.md').read()
print(len(txt))
# find section 3
import re
idx = txt.find('§3')
print(txt[:200])
