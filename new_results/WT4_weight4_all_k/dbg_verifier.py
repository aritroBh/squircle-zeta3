
import sys, os, json, re
sys.dont_write_bytecode = True
WT4HOME = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '<project>/research/08_theory/A3h4_frame_ds')
sys.path.insert(0, '<project>/research/08_theory/A3c_claude_elementary/Q')
from frame_ds import *
import cmzv as cm
import w4_symbol as ws
from qsqrt2 import QS2 as Q
print('w4_symbol file:', ws.__file__)
print('cm.A type:', type(cm.A), cm.A)
print('qsqrt2 A:', Q(1,1), cm.A == Q(1,1), hash(cm.A)==hash(Q(1,1)))
from fractions import Fraction as F
def parse_qs2(s):
    s=s.strip()
    m=re.fullmatch(r"\((.+?)([+-])(\d+(?:/\d+)?)\*sqrt2\)",s)
    if m:
        p=F(m.group(1)); q=F(m.group(3))
        if m.group(2)=='-': q=-q
        return Q(p,q)
    m=re.fullmatch(r"(-?\d+(?:/\d+)?)\*sqrt2",s)
    if m: return Q(0,F(m.group(1)))
    return Q(F(s))
def parse_word(s):
    s=s.strip(); assert s.startswith('(') and s.endswith(')')
    body=s[1:-1].strip()
    if not body: return ()
    parts,depth,cur=[],0,''
    for ch in body:
        if ch=='(': depth+=1; cur+=ch
        elif ch==')': depth-=1; cur+=ch
        elif ch==',' and depth==0: parts.append(cur); cur=''
        else: cur+=ch
    parts.append(cur)
    return tuple(parse_qs2(x) for x in parts if x.strip())
TGT='D2'
rows=[]
with open(os.path.join(WT4HOME, TGT.lower()+'.tsv')) as f:
    for line in f:
        if line.startswith('#') or line.startswith('coefficient_p'): continue
        p,q,k,w=line.rstrip('\n').split('\t')
        rows.append((parse_word(w),int(k),Q(F(p),F(q))))
tsv={}
for w,k,c in rows: tsv[(w,k)]=tsv.get((w,k),Q(0))+c
deltas=json.load(open(os.path.join(WT4HOME,'deltas_deg3.json')))
POLE={'0':Q(0),'a':cm.A,'-a':-cm.A,'b':cm.B,'-b':-cm.B}
out={}
for (i,j) in ((2,1),(3,0)):
    for name,(p,q) in deltas[f'{i},{j}'][TGT]:
        z=POLE[name]; coefficient=Q(F(p),F(q))
        for (w,k),d in ws.log_product(*((ws.LT,)*i+(ws.MT,)*j)).items():
            out[((z,)+w,k)]=out.get(((z,)+w,k),Q(0))+coefficient*d
rederived={kk:v for kk,v in out.items() if v}
print('lens', len(tsv), len(rederived))
ks = set(tsv)|set(rederived)
bad=[k for k in ks if tsv.get(k,Q(0))!=rederived.get(k,Q(0))]
print('mismatching keys:', len(bad))
for k in bad[:8]:
    print('  ', k, 'tsv=', tsv.get(k), 'rederived=', rederived.get(k))

print('bool(Q(0)):', bool(Q(0)))
zz = [k for k,v in out.items() if v==Q(0)]
print('zero-value entries in out:', len(zz))
