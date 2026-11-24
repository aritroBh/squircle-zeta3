p = 'research/08_theory/A3g_denominator/glanois_1411.4947.txt'
lines = open(p).read().split('\n')
import re
for i,l in enumerate(lines):
    if re.search(r'1 ?- ?3t|3t|Hilbert|dimension', l, re.I):
        print(i, l)
