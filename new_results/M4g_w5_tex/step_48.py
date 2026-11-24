import json, re
from fractions import Fraction as F
import mpmath as mp
mp.dps=470
pl=json.load(open('atoms_planner.json'))
r2=mp.sqrt(2)
def li(k,x):
    tot,xn,n=mp.mpf(0),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn/mp.mpf(n)**k
        if abs(xn)<mp.mpf(10)**(-455): return tot
        n+=1
def field(s):
    s=s.strip().strip('()').replace('sqrt2','s')
    v=mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+',s):
        if 's' in t:
            c=t.replace('*s','').replace('s','')
            c=c if c not in ('','+','-') else c+'1'
            v+=mp.mpf(F(c).numerator)/F(c).denominator*r2
        else:
            v+=mp.mpf(F(t).numerator)/F(t).denominator
    return v
def atom(a):
    base={'lam':mp.log(2),'L':mp.log(1+r2),'pi2':mp.pi**2,'z3':mp.zeta(3),'z5':mp.zeta(5)}
    if a in base: return base[a]
    if a=='Li2(b)': return li(2,r2-1)
    m=re.fullmatch(r'Li([2-5])\((.*)\)',a)
    return li(int(m.group(1)),field(m.group(2)))
maxd=mp.mpf(0)
for f,vs in pl.items():
    if f in ('Lchi8_5','N5'): continue
    v=atom(f)
    d=abs(v-mp.mpf(vs))
    if d>maxd: maxd=d; worst=f
    if d>mp.mpf('1e-430'):
        print('ATOM DIFF', f, mp.nstr(d,3))
print('max atom diff:', mp.nstr(maxd,3), 'at', worst)
