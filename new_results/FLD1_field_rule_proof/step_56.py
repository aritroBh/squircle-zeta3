WD='<project>/research/08_theory/FLD1_field_rule_proof'
code = r'''"""step29b_g2compare.py -- compare my k=4 with G2's B4, evaluating G2's He = c^e(H_w + (-1)^{e+nz}H_{-w}) correctly."""
import sys, json, re, pickle
from fractions import Fraction as Fr
sys.path.insert(0,'.')
from hpl import coeffs_exact
G2='<project>/research/08_theory/G2_M4alpha/'
def nz(w): return sum(1 for a in w if a!=0)
def neg(w): return tuple(-a for a in w)
def load_g2(path):
    out=[]
    for nm,v in json.load(open(path)):
        v=Fr(v)
        if nm=='1': out.append(('const',0,(),v)); continue
        if nm=='c': out.append(('const',1,(),v)); continue
        m=re.match(r'c\^(-?\d+) He?\[(.*)\]', nm)
        out.append(('He',int(m.group(1)),tuple(int(t) for t in m.group(2).split(',')),v))
    return out
def series_g2(terms,lo,hi):
    out={}
    for (typ,e,w,q) in terms:
        if typ=='const':
            if lo<=e<=hi: out[e]=out.get(e,Fr(0))+q
            continue
        sgn = 1 if (e+nz(w))%2==0 else -1
        for (ww,qq) in ((w,q),(neg(w),q*sgn)):
            if len(ww)==0: continue
            ce=coeffs_exact(ww,hi-e+2)
            for j in range(lo,hi+1):
                idx=j-e
                if 1<=idx<len(ce): out[j]=out.get(j,Fr(0))+qq*ce[idx]
    return out
def series_mine(terms,lo,hi):
    out={}
    for (e,w,q) in terms:
        if len(w)==0:
            if lo<=e<=hi: out[e]=out.get(e,Fr(0))+q
            continue
        ce=coeffs_exact(w,hi-e+2)
        for j in range(lo,hi+1):
            idx=j-e
            if 1<=idx<len(ce): out[j]=out.get(j,Fr(0))+q*ce[idx]
    return out
store = pickle.load(open('hpl_constructed_K5.pkl','rb'))
for (g2file, key, label) in [('hpl_fit_k4_(0,0,0).json',(4,(0,0,0,0)),'R4 (rational part)'),
                             ('hpl_fit_k4_(1,0,0).json',(4,(1,0,0,0)),'Psi2 (z2 part)')]:
    g2 = load_g2(G2+g2file); mine = store[key]
    sg=series_g2(g2,-4,80); sm=series_mine(mine,-4,80)
    diff={j:sg.get(j,Fr(0))-sm.get(j,Fr(0)) for j in set(sg)|set(sm)}
    diff={j:v for j,v in diff.items() if v!=0}
    print(f'{label}: G2 {len(g2)} terms vs mine {len(mine)} terms; Taylor c^-4..c^80 agree: {len(diff)==0}', ('' if not diff else str(list(diff.items())[:6])))
'''
open(WD+'/step_29b_g2compare.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_29b_g2compare.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-2000:] if r.returncode else '')