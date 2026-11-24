base = '<project>'
import re, os
d = base + '/research/08_theory/ZLAD1_zeta_ladder'
for i in list(range(4, 52)):
    f = f'{d}/step_{i:02d}.py'
    if os.path.exists(f):
        with open(f) as fh:
            head = fh.read(400)
        m = re.search(r'"""(.*?)"""', head, re.S)
        first = (m.group(1).strip().split('\n')[0] if m else head.split('\n')[0])
        print(f'step_{i:02d}: {first[:160]}')
