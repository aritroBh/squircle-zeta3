d='research/08_theory/D13h_final_cell'
print(open(d+'/step_19.py').read()[:2000])
import re
src = open(d+'/step_38.py').read()
m = re.search(r'def certified_y38.*?(?=\ndef |\nif __|\Z)', src, re.S)
print(m.group(0)[:1500] if m else "not found")
