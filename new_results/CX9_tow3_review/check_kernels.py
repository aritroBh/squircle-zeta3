"""Independent Fraction-only comparison of complete saved partial fractions."""
from fractions import Fraction as F
from pathlib import Path
from math import comb
import pickle,json
ROOT=Path(__file__).resolve().parent.parent
KEYS=((4,0),(3,1),(2,2))
H={k:pickle.load(open(ROOT/'WTS_weight_k_structure'/'theory'/f'kern_{k}.pkl','rb')) for k in range(5,10)}
J=json.loads((ROOT/'WT5_weight5_all_k'/'deltas_deg4.json').read_text())
POLES={'0':(F(0),F(0)),'a':(F(1),F(1)),'-a':(F(-1),F(-1)),'b':(F(-1),F(1)),'-b':(F(1),F(-1))}
def flat(pf):
    out={('real',pole,c):F(v) for pole,cs in pf['real'].items() for c,v in enumerate(cs)}
    out.update({(part,n,c):F(v) for part in ('ipos','ineg') for n,cs in enumerate(pf[part]) for c,v in enumerate(cs)})
    out.update({('const',0,c):F(v) for c,v in enumerate(pf['const'])})
    return out
def residual(pfs,weights):
    terms=list(map(flat,pfs)); names=set().union(*terms)
    return {n:sum(w*t.get(n,0) for w,t in zip(weights,terms)) for n in names}
nreal=0
for key in KEYS:
    pfs=[H[k][key] for k in range(5,10)]
    assert not any(residual(pfs,[1,-4,6,-4,1]).values()),key
    nreal+=4*len(set().union(*(pf['real'] for pf in pfs)))
    for k,pf in zip(range(5,10),pfs):
        assert all(F(v)==0 for part in ('ipos','ineg') for cs in pf[part] for v in cs),(key,k)
        assert all(F(v)==0 for v in pf['const']),(key,k)
    print('PASS',key,'K5=4K6-6K7+4K8-K9 as complete PF; +/-i coefficients and constant zero at k=5..9',flush=True)
    for order in (1,2,3):
        weights=[(-1)**(order-j)*comb(order,j) for j in range(order+1)]
        target={'real':{POLES[n]:tuple(map(F,cs))+(F(0),F(0)) for n,cs in J[f'{key[0]},{key[1]}'][f'D{order}']},'ipos':[],'ineg':[],'const':(0,0,0,0)}
        assert not any(residual(pfs[:order+1]+[target],weights+[-1]).values()),(key,order)
        print('PASS',key,f'Delta^{order}K5 equals WT5 D{order} (integral E{order})',flush=True)
print(f'PASS {nreal} real-residue QSI scalar components; 9 complete difference-kernel comparisons',flush=True)
