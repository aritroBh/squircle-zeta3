
import json, re
from fractions import Fraction
import mpmath as mp
mp.dps=460
print('sanity: log2 =', mp.nstr(mp.log(2),45))
print('sanity: sqrt2 =', mp.nstr(mp.sqrt(2),45))
sq2=mp.sqrt(2)
terms=json.load(open('../M4e_w5_lchi/planner/R5_nice_reduced.json'))['terms']
# ---- planner-style ----
EPS=mp.mpf(10)**(-445)
def liP(k,x):
    tot,xn,n=mp.mpf(0),mp.mpf(1),1
    while True:
        xn*=x; tot+=xn/mp.mpf(n)**k
        if abs(xn)<EPS: return tot
        n+=1
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
# ---- check_supplement-style ----
def li_series(k,x):
    N=int(450/float(-mp.log10(abs(x))))+60
    s,xp=mp.mpf(0),x
    for j in range(1,N+1):
        s+=xp/mp.mpf(j)**k; xp*=x
    return s
def parse_field_str(s):
    s=s.replace('(','').replace(')','').strip()
    p,q=Fraction(0),Fraction(0)
    for t in re.findall(r'[+-]?[^+-]+',s):
        t=t.strip()
        if t.endswith('*sqrt2'): q+=Fraction(t[:-6])
        else: p+=Fraction(t)
    return (p,q)
def mpq(fr): return mp.mpf(fr.numerator)/fr.denominator
def fval(pq): return mpq(pq[0])+mpq(pq[1])*sq2
def split_key(key):
    parts,depth,cur=[],0,''
    for ch in key:
        if ch=='*' and depth==0: parts.append(cur); cur=''; continue
        if ch=='(': depth+=1
        elif ch==')': depth-=1
        cur+=ch
    parts.append(cur)
    return parts
def atomC(key):
    v=None
    for f in split_key(key):
        m=re.fullmatch(r'Li(\d)\((.*)\)',f)
        base={'lam':mp.log(2),'L':mp.log(1+sq2),'pi2':mp.pi**2,'z3':mp.zeta(3),'z5':mp.zeta(5)}
        if f in base: fv=base[f]
        elif m:
            k=int(m.group(1)); a=m.group(2)
            x=sq2-1 if a=='b' else fval(parse_field_str(a))
            fv=li_series(k,x)
        else: raise ValueError(f)
        v=fv if v is None else v*fv
    return v
bad=0
for nm in terms:
    if nm=='Lchi8_5': continue
    d=atomP(nm)-atomC(nm)
    if abs(d)>mp.mpf('1e-400'):
        print('DIFF', nm, mp.nstr(d,3)); bad+=1
print('atoms differing:', bad)
