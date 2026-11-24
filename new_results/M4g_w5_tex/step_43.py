import json, re
from fractions import Fraction
import mpmath as mp
mp.dps=450
sq2=mp.sqrt(2)
terms=json.load(open('../M4e_w5_lchi/planner/R5_nice_reduced.json'))['terms']
n5=json.load(open('../M4e_w5_lchi/cobracket.json'))['correction']
EPS=mp.mpf(10)**(-440)
def liP(k,x):
    tot,xn,n=mp.mpf(0),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn/mp.mpf(n)**k
        if abs(xn)<EPS: return tot
        n+=1
def li32P(x,y):
    tot,H,xn,ym,n=mp.mpf(0),mp.mpf(0),mp.mpf(1),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn*H/mp.mpf(n)**3
        if n>50 and abs(xn)<EPS: return tot
        ym*=y; H+=ym/mp.mpf(n)**2; n+=1
def field(s):
    s=s.strip().strip('()').replace('sqrt2','s')
    v=mp.mpf(0)
    for t in re.findall(r'[+-]?[^+-]+',s):
        if 's' in t:
            c=t.replace('*s','').replace('s','')
            c=c if c not in ('','+','-') else c+'1'
            v+=mp.mpf(Fraction(c).numerator)/Fraction(c).denominator*sq2
        else:
            v+=mp.mpf(Fraction(t).numerator)/Fraction(t).denominator
    return v
def splitk(name):
    out,depth,cur=[],0,''
    for ch in name:
        depth+=(ch=='(')-(ch==')')
        if ch=='*' and depth==0: out.append(cur); cur=''
        else: cur+=ch
    return out+[cur]
def atomP(a):
    base={'lam':mp.log(2),'L':mp.log(1+sq2),'pi2':mp.pi**2,'z3':mp.zeta(3),'z5':mp.zeta(5)}
    if a in base: return base[a]
    if a=='Li2(b)': return liP(2,sq2-1)
    m=re.fullmatch(r'Li([2-5])\((.*)\)',a)
    return liP(int(m.group(1)),field(m.group(2)))
def mpq(fr): return mp.mpf(fr.numerator)/fr.denominator
h=[mp.zeta(5,mp.mpf(k)/8) for k in (1,3,5,7)]
Lchi=(h[0]+h[3]-h[1]-h[2])/mp.mpf(8)**5
W5c=mp.mpf(json.load(open('../M4e_w5_lchi/W5_5400.json'))['value'])
tot=mp.mpf(0)
for nm,(a,q) in terms.items():
    if nm=='Lchi8_5': v=Lchi
    else:
        v=mp.mpf(1)
        for f in splitk(nm): v*=atomP(f)
    tot+=(mpq(Fraction(a))+mpq(Fraction(q))*sq2)*v
print('classical total:',mp.nstr(tot,25))
b=sq2-1
N5v=(mpq(Fraction(779))+mpq(Fraction(-695,56))*sq2)*li32P(b,1)-672*li32P(1/sq2,1) \
    -(mpq(Fraction(374))+mpq(Fraction(-76))*sq2)*li32P(-b,1) \
    -(mpq(Fraction(405))+mpq(Fraction(-4951,56))*sq2)*li32P(b,-1)
print('N5:',mp.nstr(N5v,25))
print('total-W5:',mp.nstr(tot+N5v-W5c,4))
# now N5 from cobracket.json directly
N5b=mp.mpf(0)
for e in n5:
    c=mpq(Fraction(e['coefficient'][0]))+mpq(Fraction(e['coefficient'][1]))*sq2
    x=mpq(Fraction(e['x'][0]))+mpq(Fraction(e['x'][1]))*sq2
    y=mpq(Fraction(e['y'][0]))+mpq(Fraction(e['y'][1]))*sq2
    N5b+=c*li32P(x,y)
print('N5 from json:',mp.nstr(N5b,25))
print('N5 json vs hard:',mp.nstr(N5b-N5v,4))
