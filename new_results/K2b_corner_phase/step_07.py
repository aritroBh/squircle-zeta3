
wd = 'research/08_theory/K2b_corner_phase/'
utils = '''import re, mpmath
mpmath.mp.dps = 80
BASE = 'research/08_theory/K2_corner_decay/out/'
FILES = {'pi_3':'mk_pi_3.tsv','pi_2':'mk_pi_2.tsv','2pi_3':'mk_2pi_3.tsv','3pi_4':'mk_3pi_4.tsv','5pi_6':'mk_5pi_6.tsv'}
ALPHA = {'pi_3': mpmath.pi/3, 'pi_2': mpmath.pi/2, '2pi_3': 2*mpmath.pi/3, '3pi_4': 3*mpmath.pi/4, '5pi_6': 5*mpmath.pi/6}

def parse_complex(s):
    s = s.strip().strip('()')
    idx = None
    for i in range(len(s)-1, 0, -1):
        if s[i] in '+-' and s[i-1] not in 'eE':
            idx = i; break
    return mpmath.mpf(s[:idx]), mpmath.mpf(s[idx:].replace('j','').replace('i',''))

def load(key):
    d = {}
    with open(BASE+FILES[key]) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split('\\t')
            if len(parts) < 2: parts = line.split()
            r,i = parse_complex(parts[1])
            d[int(parts[0])] = mpmath.mpc(r,i)
    return d
'''
open(wd+'k2utils.py','w').write(utils)
# inspect raw lines
base = 'research/08_theory/K2_corner_decay/out/'
with open(base+'mk_pi_2.tsv') as f:
    lines = f.readlines()
print(len(lines))
for ln in lines[:2] + lines[100:102] + lines[299:301]:
    print(repr(ln))
