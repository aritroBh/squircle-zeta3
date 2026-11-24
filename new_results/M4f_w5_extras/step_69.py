import subprocess, re
WORK='<project>/research/08_theory/M4f_w5_extras'
# W4pp piece vs M4a's W4pp_800.txt, and per-piece totals
src = r'''
import json, re, sys
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
HERE = Path("<project>/research/08_theory/M4f_w5_extras")
sys.path.insert(0, str(HERE))
mp.mp.dps = 880
import importlib.util
spec = importlib.util.spec_from_file_location("cm", HERE/"check_m4.py")
# instead of importing (it runs), copy needed funcs inline:
r2, lam, L = mp.sqrt(2), mp.log(2), mp.log(1+mp.sqrt(2))
def field(s):
    s = s.strip()
    if s == 'b': return r2 - 1
    if s == '-b': return 1 - r2
    if s == '1/sqrt2': return 1/r2
    s = s.strip('()').replace('sqrt2','s')
    v = mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+', s):
        if 's' in t:
            c = t.replace('*s','').replace('s','')
            c = c if c not in ('','+','-') else c+'1'
            v += mp.mpf(F(c).numerator)/F(c).denominator*r2
        else:
            v += mp.mpf(F(t).numerator)/F(t).denominator
    return v
cache={}
def atom(a):
    if a in cache: return cache[a]
    if a in ('lam','log2'): v = lam
    elif a == 'L': v = L
    elif a in ('pi2','pi^2'): v = mp.pi**2
    elif a == 'z3': v = mp.zeta(3)
    else:
        m = re.fullmatch(r'Li(\d)\((.*)\)', a); v = mp.polylog(int(m.group(1)), field(m.group(2)))
    cache[a]=v; return v
def split(nm):
    out, depth, cur = [], 0, ''
    for ch in nm:
        if ch=='(' : depth+=1
        if ch==')' : depth-=1
        if ch=='*' and depth==0: out.append(cur); cur=''
        else: cur+=ch
    out.append(cur); return out
def cf(p,q): return mp.mpf(F(p).numerator)/F(p).denominator + r2*mp.mpf(F(q).numerator)/F(q).denominator
def prodval(key):
    key = key.replace('log2^3','log2*log2*log2').replace('log2^2','log2*log2').replace('L^3','L*L*L').replace('L^2','L*L')
    v = mp.mpf(1)
    for f in split(key): v *= atom(f)
    return v
d = json.loads((HERE/'m4_terms.json').read_text())
W4 = mp.mpf(0)
for k,(p,q) in d['W4pp'].items(): W4 += cf(p,q)*prodval(k)
W34 = mp.mpf(0)
for k,(p,q) in d['W34'].items(): W34 += cf(p,q)*prodval(k)
T = mp.mpf(0)
for k,(p,q) in d['T'].items(): T += cf(p,q)*prodval(k)
w4ref = mp.mpf(re.sub(r'\s','',(HERE.parent/'M4a_decomposition'/'W4pp_800.txt').read_text()).split('=')[-1])
print('W4pp - W4pp_800 =', mp.nstr(W4-w4ref, 8))
print('W34 =', mp.nstr(W34, 30))
print('T =', mp.nstr(T, 30))
print('W4pp =', mp.nstr(W4, 30))
'''
open(f'{WORK}/dbg_m4.py','w').write(src)
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','dbg_m4.py'],
                   cwd=WORK, capture_output=True, text=True, timeout=600)
print(r.stdout, r.stderr[-300:])
