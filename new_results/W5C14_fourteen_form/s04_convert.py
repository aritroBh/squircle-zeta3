"""Step 3 (exact): convert R14 - R5 into exact level-8 H-word vectors (rational and sqrt2 parts).
Uses the project's audited conversion code (cmzv.forms_path / atom), as W5C1's w5_build_D.py.
New atoms (not cross-ratios of S):
  Li5(b^4)        : t=(v/a)^2,               v: 0 -> b,   letters 0, +-a
  Li5((4-3r2)/8)  : t=s/(s-1), s=(v/a)^2,    v: 0 -> b    ((4-3r2)/8 = b^4/(b^4-1))
  Li5(1-b^4)      : t=1-s,     s=(v/a)^2,    v: a -> b
Li5(x) = int_0^x (dt/t)^4 o dt/(1-t); forms listed outer -> inner (cmzv convention, checked numerically in s05).
Output: delta.pkl {rat:{word:F}, sq:{word:F}, atoms:{name:vec}}"""
import sys, pickle, time
sys.dont_write_bytecode=True
TH='<project>/research/08_theory'
sys.path.insert(0, TH+'/A3h4_frame_ds'); sys.path.insert(0, TH+'/A3c_claude_elementary/Q')
from fractions import Fraction as F
from collections import Counter
import cmzv as cm
from cmzv import Q, POINTS, dlog, forms_path, add, mul
from common import *
A,B=cm.A,cm.B
t0=time.monotonic()
def comb(*pc):
    """linear combination of dlog(x-p): pc = (p,c),...  -> cmzv form"""
    out=[]
    for p,c in pc:
        for e,s in dlog(p): out.append((e,s*c))
    return tuple(out)
def path_neg(forms,c,e): return {w:-x for w,x in forms_path(forms,POINTS[c],POINTS[e]).items()}
NEW={}
z0=Q(0)
f_v=comb((z0,2)); f_pm=comb((A,1),(-A,1))
NEW[(F(17),F(-12))]       = forms_path([f_v]*4+[comb((A,-1),(-A,-1))], POINTS[z0], POINTS[B])
NEW[(F(1,2),F(-3,8))]     = forms_path([comb((z0,2),(A,-1),(-A,-1))]*4+[f_pm], POINTS[z0], POINTS[B])
NEW[(F(-16),F(12))]       = forms_path([f_pm]*4+[comb((z0,-2))], POINTS[A], POINTS[B])
def atom_vec(k):
    t=k[0]
    if t=='lam': return cm.atom(('log',Q(2)))
    if t=='L': return cm.atom(('log',A))
    if t=='pi2': return cm.atom(('pi2',))
    if t=='z3': return cm.atom(('z3',))
    if t=='z5': return {(-1,-1,-1,-1,0):F(1)}
    _,n,(p,q)=k
    if n==5 and (p,q) in NEW: return NEW[(p,q)]
    return cm.atom(('li',n,Q(p,q)))
LCHI={(-1,-1,-1,-1,1):F(1,4),(-1,-1,-1,-1,3):F(-1,4),(-1,-1,-1,-1,5):F(-1,4),(-1,-1,-1,-1,7):F(1,4)}
def form_vec(form):
    rat,sq=Counter(),Counter()
    for key,(cp,cq) in terms_keys(form).items():
        if key==(('Lchi',),):
            for w,a in LCHI.items(): rat[w]+=2*cq*a; sq[w]+=cp*a
            continue
        v=None
        for a in key:
            av=atom_vec(a); v=dict(av) if v is None else mul(v,av)
        for w,a in v.items(): rat[w]+=cp*a; sq[w]+=cq*a
    return ({w:a for w,a in rat.items() if a},{w:a for w,a in sq.items() if a})
r14=form_vec(load(R14_PATH)); r5=form_vec(load(R5_PATH))
d=[]
for j in range(2):
    c=Counter(r14[j])
    for w,a in r5[j].items(): c[w]-=a
    d.append({w:a for w,a in c.items() if a})
for name,v in [('R14',r14),('R5',r5),('Delta',d)]:
    print(f'{name}: support rat/sq = {len(v[0])}/{len(v[1])}')
bad=[w for w in list(d[0])+list(d[1]) if w[0]==0 or w[-1]==-1]
print('non-admissible words in Delta:',len(bad))
H=lambda v: max((max(abs(a.numerator),a.denominator) for a in v.values()),default=0)
print('Delta heights rat/sq:',H(d[0]),H(d[1]))
# which Li5 atoms appear in Delta with nonzero coefficient
k14=terms_keys(load(R14_PATH)); k5=terms_keys(load(R5_PATH))
diff={k:(k14.get(k,(0,0))[0]-k5.get(k,(0,0))[0], k14.get(k,(0,0))[1]-k5.get(k,(0,0))[1]) for k in set(k14)|set(k5)}
diff={k:v for k,v in diff.items() if v!=(0,0)}
print('classical terms in R14-R5:',len(diff))
for k,v in sorted(diff.items(),key=str):
    if any(a[0]=='Li' and a[1]==5 for a in k) or k in ((('z5',),),(('Lchi',),)): print('  ',k,v)
pickle.dump({'rat':d[0],'sq':d[1],'new_atoms':NEW,'r14':r14,'r5':r5,'diff_terms':diff},open('delta.pkl','wb'))
print(f'seconds={time.monotonic()-t0:.1f}')
