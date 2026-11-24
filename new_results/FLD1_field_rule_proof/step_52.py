WD='<project>/research/08_theory/FLD1_field_rule_proof'
code = r'''"""step29_g2compare.py -- (a) is the UNSYMMETRISED construction already even? (b) compare k=4 with G2's B4."""
import sys, json, re, pickle
from fractions import Fraction as Fr
sys.path.insert(0,'.')
import zsum
from zsum import K, add_into
from hpl import coeffs_exact
# rebuild unsymmetrised terms via construct()
src = open("step_25_construct.py").read()
ns={}
exec(src[:src.index("if __name__")], ns)
Zin = pickle.load(open('mu_zsum_K5.pkl','rb'))
construct=ns['construct']; symmetrise=ns['symmetrise']; odd_form=ns['odd_form']; even_form=ns['even_form']
def oddlen(terms):
    return len(odd_form(terms))
for k in (2,3,4):
    for zkey,E in sorted(Zin[k].items()):
        terms = construct(E)
        print(f'k={k} z={zkey}: unsym {len(terms)} terms, odd-form words of UNSYM = {oddlen(terms)}')
# (b) compare k=4 with G2
G2='<project>/research/08_theory/G2_M4alpha/'
def load_fit(path):
    out=[]
    for nm,v in json.load(open(path)):
        if nm=='1': out.append((0,(),Fr(v))); continue
        if nm=='c': out.append((1,(),Fr(v))); continue
        m=re.match(r'c\^(-?\d+) H\[(.*)\]', nm); out.append((int(m.group(1)),tuple(int(t) for t in m.group(2).split(',')),Fr(v)))
    return out
def series(terms,lo,hi):
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
# G2 R4 (rational part) vs my k=4 rational part (symmetrised+corrected)
g2R4 = load_fit(G2+'hpl_fit_k4_(0,0,0).json')
mine = store[(4,(0,0,0,0))]
sg=series(g2R4,-4,60); sm=series(mine,-4,60)
diff={j:sg.get(j,Fr(0))-sm.get(j,Fr(0)) for j in set(sg)|set(sm)}
diff={j:v for j,v in diff.items() if v!=0}
print()
print(f'G2 R4 ({len(g2R4)} terms) vs my k4 rational part ({len(mine)} terms): Taylor c^-4..c^60 agree: {len(diff)==0}', ('' if not diff else str(list(diff.items())[:5])))
# G2 psi2 (z2 part) vs my k=4 z2 part
g2psi2 = load_fit(G2+'hpl_fit_k4_(1,0,0).json')
mine2 = store[(4,(1,0,0,0))]
sg=series(g2psi2,-4,60); sm=series(mine2,-4,60)
diff={j:sg.get(j,Fr(0))-sm.get(j,Fr(0)) for j in set(sg)|set(sm)}
diff={j:v for j,v in diff.items() if v!=0}
print(f'G2 Psi2 ({len(g2psi2)} terms) vs my k4 z2 part ({len(mine2)} terms): agree: {len(diff)==0}', ('' if not diff else str(list(diff.items())[:5])))
'''
open(WD+'/step_29_g2compare.py','w').write(code)
import subprocess
VENV='<project>/.venv/bin/python'
r = subprocess.run(['nice',VENV,'step_29_g2compare.py'],capture_output=True,text=True,cwd=WD)
print(r.stdout); print(r.stderr[-3000:] if r.returncode else '')