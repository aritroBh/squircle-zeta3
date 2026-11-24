import json, re
from fractions import Fraction as F
import mpmath as mp
mp.dps=470
r2=mp.sqrt(2); lam=mp.log(2); L=mp.log(1+r2)
EPS=mp.mpf(10)**(-460)
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
def li(k,x):
    if x==-1: return -(1-mp.mpf(2)**(1-k))*mp.zeta(k)
    tot,xn,n=mp.mpf(0),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn/mp.mpf(n)**k
        if abs(xn)<EPS: return tot
        n+=1
def li32(x,y):
    tot,H,xn,ym,n=mp.mpf(0),mp.mpf(0),mp.mpf(1),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn*H/mp.mpf(n)**3
        if n>50 and abs(xn)<EPS: return tot
        ym*=y; H+=ym/mp.mpf(n)**2; n+=1
def atom(a):
    base={'lam':lam,'L':L,'pi2':mp.pi**2,'z3':mp.zeta(3),'z5':mp.zeta(5)}
    if a in base: return base[a]
    if a=='Li2(b)': return li(2,r2-1)
    m=re.fullmatch(r'Li([2-5])\((.*)\)',a); assert m,a
    return li(int(m.group(1)),field(m.group(2)))
def split(name):
    out,depth,cur=[],0,''
    for ch in name:
        depth+=(ch=='(')-(ch==')')
        if ch=='*' and depth==0: out.append(cur); cur=''
        else: cur+=ch
    return out+[cur]
def num(p,q): return mp.mpf(p.numerator)/p.denominator+r2*mp.mpf(q.numerator)/q.denominator
terms=json.load(open('../M4e_w5_lchi/planner/R5_nice_reduced.json'))['terms']
b=r2-1
N5=(num(F(779),F(-695,56))*li32(b,1)-672*li32(1/r2,1)
    -num(F(374),F(-76))*li32(-b,1)-num(F(405),F(-4951,56))*li32(b,-1))
h=[mp.zeta(5,mp.mpf(k)/8) for k in (1,3,5,7)]
Lchi=(h[0]+h[3]-h[1]-h[2])/mp.mpf(8)**5
tot=N5
for nm,(a,q) in terms.items():
    if nm=='Lchi8_5': v=Lchi
    else:
        v=mp.mpf(1)
        for f in split(nm): v*=atom(f)
    tot+=num(F(a),F(q))*v
W5c=mp.mpf(json.load(open('../M4e_w5_lchi/W5_5400.json'))['value'])
print('planner-style eval of REDUCED json - W5 =', mp.nstr(tot-W5c,6))
print('N5 =', mp.nstr(N5,30))
