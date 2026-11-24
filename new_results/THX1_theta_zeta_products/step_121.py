
src = open('step_115.py').read()
import re
# find X8 and LI definitions in step_115
for line in src.split('\n'):
    if any(k in line for k in ['X8c', 'X8d', 'LI9', 'X8a', 'X8b', 'def li_half']):
        print(line)
