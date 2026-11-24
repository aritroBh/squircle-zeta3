import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import mpmath as mp
import numeric_cmzv as NC
from hexc2_common import hpl_to_cmzv, li_atom
from fractions import Fraction as Fr
sys.path.insert(0,'<project>/research/08_theory/COR4_M3alpha_classical')
from hpl_series import H_series

DPS = 70
NC.set_dps(DPS)
import json, re
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle'
def load_fit(name):
    D = json.load(open(G1+'/'+name))
    out=[]
    for key,val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key=="1": e,w = 0,()
        elif key=="c": e,w = 1,()
        else:
            e=int(m.group(1)); w=tuple(int(z) for z in m.group(2).split(',')) if m.group(2).strip() else ()
        out.append((e,w))
    return out
WORDS = sorted({w for name in ['hpl_fit_k3_(0,0,0).json','hpl_fit_k3_(1,0,0).json'] for e,w in load_fit(name) if w}, key=lambda x:(len(x),x))
print('distinct HPL words:', len(WORDS))
def cmzv_combo(d):
    return sum(mp.mpc(int(c.numerator), int(c.denominator)).real/ c.denominator * NC.cmzv_value(w) for w,c in d.items())
worst = 0; fails = []
for w in WORDS:
    conv = hpl_to_cmzv(w)
    got = cmzv_combo(conv)
    expect = H_series(w, mp.mpf(1)/2)
    err = abs(got-expect)
    rel = err/max(1,abs(expect))
    worst = max(worst, float(rel))
    if rel > 1e-62: fails.append((w, float(rel)))
print('worst relative err over %d words: %.2e'%(len(WORDS), worst))
print('fails:', fails[:10])
