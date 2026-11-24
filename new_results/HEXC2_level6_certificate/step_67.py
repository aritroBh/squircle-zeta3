import sys, time, json, re
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import mpmath as mp
import fast_cmzv as FC
from hexc2_common import hpl_to_cmzv, li_atom
from hpl_series import H_series
from fractions import Fraction as Fr
FC.setup(62); FC.finish_setup()
mp.mp.dps = 65
t0 = time.monotonic()
def combo_fast(d):
    tot = mp.mpc(0)
    for w,c in d.items():
        v = FC.cmzv_value(w)
        tot += mp.mpf(c.numerator)/c.denominator * mp.mpc(v.real.mid(), v.imag.mid())
    return tot
G1 = '<project>/research/10_frontier/G1_zeta3_any_angle'
def load_fit(name):
    D = json.load(open(G1+'/'+name)); out=[]
    for key,val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key=="1": e,w = 0,()
        elif key=="c": e,w = 1,()
        else: e=int(m.group(1)); w=tuple(int(z) for z in m.group(2).split(',')) if m.group(2).strip() else ()
        out.append((e,w))
    return out
WORDS = sorted({w for name in ['hpl_fit_k3_(0,0,0).json','hpl_fit_k3_(1,0,0).json'] for e,w in load_fit(name) if w}, key=lambda x:(len(x),x))
worst=0
for w in WORDS:
    got = combo_fast(hpl_to_cmzv(w))
    err = abs(got - H_series(w, mp.mpf(1)/2))
    worst = max(worst, float(err))
t1 = time.monotonic()
print(f'38 HPL conversions: {t1-t0:.1f}s worst err {worst:.2e}')
LI = [(2,'1/3'),(3,'1/3'),(3,'-1/2'),(4,'1/3'),(4,'-1/3'),(4,'2/3'),(4,'1/2'),(4,'-1/2'),(4,'-2')]
worst2=0
for n,z in LI:
    got = combo_fast(li_atom(n, Fr(z)))
    err = abs(got - mp.polylog(n, mp.mpf(z)))
    worst2 = max(worst2, float(err))
t2 = time.monotonic()
print(f'9 Li atoms: {t2-t1:.1f}s worst err {worst2:.2e}')
