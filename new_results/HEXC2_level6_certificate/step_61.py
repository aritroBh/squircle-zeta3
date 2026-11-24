import sys, time, json, re
sys.dont_write_bytecode = True
sys.path.insert(0, '<project>/research/08_theory/COR4_M3alpha_classical')
sys.path.insert(0, '<project>/research/08_theory/HEXC1_hexagon_certificate')
sys.path.insert(0, '.')
import mpmath as mp
import numeric_cmzv as NC
from hexc2_common import hpl_to_cmzv, li_atom, atom_log2, atom_log3, atom_ipi, atom_zeta2, atom_zeta3, atom_pi2, atom_pi4
from hpl_series import H_series
from fractions import Fraction as Fr
NC.set_dps(62)
t0 = time.monotonic()
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
def combo(d):
    tot = mp.mpc(0)
    for w,c in d.items():
        tot += mp.mpf(c.numerator)/c.denominator * NC.cmzv_value(w)
    return tot
worst = 0
for w in WORDS:
    got = combo(hpl_to_cmzv(w))
    err = abs(got - H_series(w, mp.mpf(1)/2))
    worst = max(worst, float(err))
t1 = time.monotonic()
print(f'38 HPL conversions at 62dps: {t1-t0:.1f}s, worst abs err {worst:.2e}')
LI = [(2,'1/3'),(3,'1/3'),(3,'-1/2'),(4,'1/3'),(4,'-1/3'),(4,'2/3'),(4,'1/2'),(4,'-1/2'),(4,'-2')]
worst2 = 0
for n,z in LI:
    got = combo(li_atom(n, Fr(z)))
    err = abs(got - mp.polylog(n, mp.mpf(z)))
    worst2 = max(worst2, float(err))
t2 = time.monotonic()
print(f'9 Li atoms: {t2-t1:.1f}s, worst abs err {worst2:.2e}')
# elementary atoms
z3v = combo(atom_zeta3()); z2v = combo(atom_zeta2())
p2v = combo(atom_pi2()); p4v = combo(atom_pi4())
print('zeta3 err %.2e, zeta2 err %.2e, pi2 err %.2e (=val -(pi)^2? check), pi4 err %.2e'%(
    abs(z3v-mp.zeta(3)), abs(z2v-mp.zeta(2)), abs(p2v+mp.pi**2), abs(p4v-mp.pi**4)))
print('note: atom_pi2 value = (i pi)^2 = -pi^2 by construction; residual above vs -pi^2')
