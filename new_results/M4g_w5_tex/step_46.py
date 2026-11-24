import json, re
from fractions import Fraction as F
import mpmath as mp
terms=json.load(open('../M4e_w5_lchi/planner/R5_nice_reduced.json'))['terms']
def run(dps):
    mp.dps=dps
    r2=mp.sqrt(2); EPS=mp.mpf(10)**(-(dps-15))
    def li(k,x):
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
                v+=mp.mpf(F(c).numerator)/F(c).denominator*r2
            else:
                v+=mp.mpf(F(t).numerator)/F(t).denominator
        return v
    def splitk(name):
        out,depth,cur=[],0,''
        for ch in name:
            depth+=(ch=='(')-(ch==')')
            if ch=='*' and depth==0: out.append(cur); cur=''
            else: cur+=ch
        return out+[cur]
    def atom(a):
        base={'lam':mp.log(2),'L':mp.log(1+r2),'pi2':mp.pi**2,'z3':mp.zeta(3),'z5':mp.zeta(5)}
        if a in base: return base[a]
        if a=='Li2(b)': return li(2,r2-1)
        m=re.fullmatch(r'Li([2-5])\((.*)\)',a)
        return li(int(m.group(1)),field(m.group(2)))
    def num(p,q): return mp.mpf(p.numerator)/p.denominator+r2*mp.mpf(q.numerator)/q.denominator
    h=[mp.zeta(5,mp.mpf(k)/8) for k in (1,3,5,7)]
    Lchi=(h[0]+h[3]-h[1]-h[2])/mp.mpf(8)**5
    groups={'li5':[],'li4':[],'w3w2':[],'pure':[]}
    for nm,(a,q) in terms.items():
        g='li5' if (nm.startswith('Li5') or nm in ('z5','Lchi8_5')) else ('li4' if 'Li4' in nm else ('pure' if set(splitk(nm))<= {'lam','L','pi2'} else 'w3w2'))
        if nm=='Lchi8_5': v=Lchi
        else:
            v=mp.mpf(1)
            for f in splitk(nm): v*=atom(f)
        groups[g].append(num(F(a),F(q))*v)
    return {g:mp.fsum(v) for g,v in groups.items()}
g1=run(450); g2=run(600)
for g in g1:
    d=g1[g]-g2[g]
    print(g, mp.nstr(g1[g],22), mp.nstr(d,3))
