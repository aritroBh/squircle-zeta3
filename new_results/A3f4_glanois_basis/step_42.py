p = 'research/08_theory/A3g_denominator/glanois_1411.4947.txt'
lines = open(p).read().split('\n')
import re
hits = [i for i,l in enumerate(lines) if '4.4' in l or 'case N = 4' in l or 'case N=4' in l]
print(hits)
# find section 4.4.2
for i,l in enumerate(lines):
    if 'N = 4' in l and 1900 < i < 2200:
        print(i, l)
