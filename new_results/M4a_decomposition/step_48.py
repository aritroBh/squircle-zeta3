import os
BASE = '<project>'
OUT = os.path.join(BASE, 'paper/squircle_note/outputs')
raw = open(os.path.join(OUT,'M3_methodA_1200.txt')).read()
print(repr(raw[:100]))
print('len', len(raw.strip()))
import mpmath as mp
mp.mp.dps = 1300
M3file = mp.mpf(raw.strip())
print(mp.nstr(M3file, 70))
# my 200-digit quad gave ...3790969819666999; closed form ...379096981967
