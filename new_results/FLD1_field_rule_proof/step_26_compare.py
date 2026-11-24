"""step26_compare.py -- compare constructed HPL forms (step_25) with the known hpl_fit forms (G1) for k=2,3.
Two forms may differ by HPL shuffle identities; we compare by evaluating the DIFFERENCE as a power series
(exact Taylor coefficients): if identically zero, the forms are equal as functions."""
import sys, json, re, pickle
from fractions import Fraction as Fr
sys.path.insert(0,'.')
from hpl import coeffs_exact
G1='<project>/research/10_frontier/G1_zeta3_any_angle/'
store = pickle.load(open('hpl_constructed_K3.pkl','rb'))
def load_fit(name):
    out=[]
    for nm,v in json.load(open(G1+name)):
        if nm=='1': out.append((0,(),Fr(v))); continue
        if nm=='c': out.append((1,(),Fr(v))); continue
        m=re.match(r'c\^(-?\d+) H\[(.*)\]', nm)
        out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), Fr(v)))
    return out
FIT = {(2,(0,0,0,0)):'hpl_fit_k2_(0,0,0).json',(3,(0,0,0,0)):'hpl_fit_k3_(0,0,0).json',
       (3,(1,0,0,0)):'hpl_fit_k3_(1,0,0).json',(3,(0,1,0,0)):'hpl_fit_k3_(0,1,0).json'}
def series_coeffs(terms, lo, hi):
    out={}
    for (e,w,q) in terms:
        if len(w)==0:
            if lo<=e<=hi: out[e]=out.get(e,Fr(0))+q
            continue
        ce = coeffs_exact(w, hi-e+2)
        for j in range(lo,hi+1):
            idx=j-e
            if 1<=idx<len(ce): out[j]=out.get(j,Fr(0))+q*ce[idx]
    return out
allok=True
for key,fn in FIT.items():
    fit = load_fit(fn)
    mine = store.get(key, [])
    cf = series_coeffs(fit,-3,40); cm = series_coeffs(mine,-3,40)
    jj = set(cf)|set(cm)
    diff = {j: cf.get(j,Fr(0))-cm.get(j,Fr(0)) for j in jj}
    diff = {j:v for j,v in diff.items() if v!=0}
    ok = len(diff)==0
    allok = allok and ok
    print(f'{key}: fit {len(fit)} terms, mine {len(mine)} terms; Taylor coeffs c^-3..c^40 agree: {ok}', ('' if ok else ('DIFF '+str(list(diff.items())[:5]))))
print('constructed k=2,3 forms == fitted B2/B3 forms (as functions):', allok)
