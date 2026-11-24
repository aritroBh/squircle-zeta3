import sys, json, pickle
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
import os
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
import mpmath as mp
mp.mp.dps = 100
from w6_lll import build, col_values
d = build(50)
nk = len(d['kv'])
vals = col_values(d, list(range(2*nk)))
logs = [mp.log10(abs(v)) if v != 0 else mp.nan for v in vals]
S = mp.fsum(l for l in logs[1:] if l == l)
print('n cols =', len(vals), ' S = sum log10|cols| =', mp.nstr(S, 8))
print('target =', mp.nstr(d['target'], 20))
print('min/max log10|kv|:', mp.nstr(min(logs[1:]),6), mp.nstr(max(logs[1:]),6))
for D in [1300, 1550, 1700]:
    print(f'det^(1/n) at {D} digits ~ 10^{mp.nstr((D+float(S))/179, 7)}')
