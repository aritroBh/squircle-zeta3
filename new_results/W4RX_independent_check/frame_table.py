"""W4RX: the 28 frame rows written out (frame, points, identity, convergence) -> logs + LaTeX rows."""
import sys, pickle
sys.dont_write_bytecode = True
from collections import OrderedDict
from w4rx_core import *
b = pickle.load(open('data/check_b.pickle', 'rb'))
NAMES = {repr(K(0, F(1, 2))): r'\tfrac1{\sqrt2}', repr(K(0, F(-1, 2))): r'-\tfrac1{\sqrt2}',
         repr(K(2, -1)): r'2-\sqrt2', repr(K(-2, 2)): r'2b', repr(K(1, -1)): r'-b',
         repr(K(1, F(-1, 2))): r'1-\tfrac1{\sqrt2}', repr(K(-1, 1)): 'b'}
PT = {repr(K(0)): '0', repr(K(1)): '1', repr(K(-1)): '-1', repr(A): 'a', repr(-A): '-a', repr(B): 'b', repr(-B): '-b', 'None': r'\infty'}
def nm(z): return NAMES.get(repr(z), repr(z))
def pt(p): return PT[repr(p)]
groups = OrderedDict()
for t in b['frame_table']:
    groups.setdefault((t['n'], t['i']), []).append(t)
lines = []
for (n, i), ts in groups.items():
    t = ts[0]
    p0, p1, pinf = t['frame']
    prod = r'\,'.join(r'\Li_{%s}(%s)' % (','.join(str(s) for s, z in f), ','.join(nm(z) for s, z in f)) for f in t['factors'])
    lifts = ', '.join('(' + ','.join(str(x) for x in s['lift']) + ')' for s in ts) if t['lift'] else '---'
    print(f'w{n} #{i}: frame (p0,p1,pinf)=({pt(p0)},{pt(p1)},{pt(pinf)}); product {prod}; '
          f'max |cumulative product| = {t["maxmod"]} = {float(t["maxmod"].mpf()):.4f}; identity words {t["ident_terms"]}; lifts {lifts}')
    print('    frame coordinates t(x):', {PT[k] if k in PT else k: v for k, v in t['ts'].items()})
    lines.append(f'{n} & {i} & $({pt(p0)},{pt(p1)},{pt(pinf)})$ & ${prod}$ & ${nm(t["maxmod"])}$ & {len(ts)} \\\\')
open('data/frame_table.tex', 'w').write('\n'.join(lines) + '\n')
