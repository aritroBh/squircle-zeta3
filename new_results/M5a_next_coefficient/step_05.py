base = '<project>'
print(open(f'{base}/research/08_theory/M5a_next_coefficient.brief').read())
print('='*80)
import os
p = f'{base}/research/08_theory/M5a_next_coefficient'
print(os.listdir(p) if os.path.exists(p) else "folder does not exist")
