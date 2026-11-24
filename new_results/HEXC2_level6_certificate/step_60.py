import sys, time, json, re
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import mpmath as mp
import numeric_cmzv as NC
from hexc2_common import hpl_to_cmzv, li_atom, target2, target3
from fractions import Fraction as Fr
sys.path.insert(0,'<project>/research/08_theory/COR4_M3alpha_classical')
from hpl_series import H_series
NC.set_dps(62)
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle'
def load_fit(name):
    D = json.load(open(G1+'/'+name))
    out=[]
    for key,val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if not m: continue
        e=int(m.group(1)); w=tuple(int(z) for z in m.group(2).split(',')) if m.group(2).strip() else ()
        out.append((e,w))
    return out
WORDS = sorted({w for name in ['hpl_fit_k3_(0,0,0).json','hpl_fit_k3_(1,0,0).json'] for e,w in load_fit(name) if w}, key=lambda x:(len(x),x))
def combo(d):
    tot = mp.mpc(0)
    for w,c in d.items():
        tot += mp.mpf(c.numerator)/c.denominator * NC.cmzv_value(w)
    return tot
t0=time.monotonic()
worst=0
for w in WORDS:
    got = combo(hpl_to_cmzv(w)); expect = H_series(w, mp.mpf(1)/2)
    rel = abs(got-expect)/max(1,abs(expect)); worst=max(worst,float(rel))
t_hpl = time.monotonic()-t0
print('38 HPL words at 62dps: %.1fs, worst rel err %.2e'%(t_hpl, worst))
t0=time.monotonic()
worst2=0
for n,zs in {2:[Fr(1,3)],3:[Fr(1,3),Fr(-1,2)],4:[Fr(1,3),Fr(-1,3),Fr(2,3),Fr(1,2),Fr(-1,2),Fr(-2)]}.items():
    for z in zs:
        got = combo(li_atom(n,z)); expect = mp.polylog(n,z)
        rel = abs(got-expect)/max(1,abs(expect)); worst2=max(worst2,float(rel))
t_li = time.monotonic()-t0
print('9 Li atoms: %.1fs, worst %.2e'%(t_li, worst2))
t0=time.monotonic()
v2 = combo(target2()); v3 = combo(target3())
print('T2+T3 values: %.1fs, residuals %s %s'%(time.monotonic()-t0, mp.nstr(v2,4), mp.nstr(v3,4)))
