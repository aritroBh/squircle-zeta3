
import os
wd = 'research/08_theory/K2b_corner_phase/'
print(open(wd+'result.md').read()[:3000])
print('====PROGRESS LOG====')
try:
    print(open(wd+'progress.log').read()[-3000:])
except Exception as e:
    print('no progress.log', e)
