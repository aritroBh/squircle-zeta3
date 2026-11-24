p = 'research/08_theory/A3g_denominator/glanois_1411.4947.txt'
with open(p) as f:
    txt = f.read()
lines = txt.split('\n')
# find where multiple zeta relative to mu_N is defined
import re
for i, l in enumerate(lines):
    if 'multiple zeta values relative' in l.lower() or 'relative to' in l.lower():
        print(i, l)
