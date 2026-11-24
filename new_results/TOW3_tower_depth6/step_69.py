import sys, os, pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
tn = f'{base}/TOWN1_tower_large_n'
wts = f'{base}/WTS_weight_k_structure/theory'
sys.path.insert(0, tn); sys.path.insert(0, wts)
os.chdir(tn)
import kernel_from_tower as KFT
from fastred import QSI, ZERO, ONE, IQ, AQ, BQ
# build tower dict {r: QSI} for key (3,0) degree n=3 at level k=6 from WTS tower_6.pkl
T6 = pickle.load(open(f'{wts}/tower_6.pkl','rb'))
rec = T6.get((3,0), {})
h = {r: QSI(*rec[r]) for r in rec if r >= 2}
print('orders:', sorted(h.keys()))
gm = KFT.gm_from_tower(h)
# Gm(-a): evaluate PF gm at v = -a. gm = [real, ipos, ineg, const]
# ipos[m-1] = e_m (order m). Gm(-a) = sum_m e_m (-a-i)^{-m} + econj_m (-a+i)^{-m} + const
import sympy as sp
s2 = sp.sqrt(2); ii = sp.I
av = 1+s2; bv = s2-1
tot = gm[3]
def qsi_sp(q): return sp.Rational(q.p0) + sp.Rational(q.p1)*s2 + sp.Rational(q.p2)*ii + sp.Rational(q.p3)*s2*ii
for idx, e in enumerate(gm[1]):
    m_ = idx  # ipos[0] is order 1 (zero), ipos[m-1] = order m
    if idx == 0: continue
    tot += qsi_sp(e) * (-av-ii)**(-idx)
for idx, e in enumerate(gm[2]):
    if idx == 0: continue
    tot += qsi_sp(e) * (-av+ii)**(-idx)
tot = sp.expand(tot)
print('Gm(-a) via kernel_from_tower:', sp.simplify(tot))
# compare: Res_{-a} K_{2,0}^{(6)} = -3 * Gm(-a) should equal recorded kern_6[(2,0)]['real'][-a]
K6 = pickle.load(open(f'{wts}/kern_6.pkl','rb'))
rec = K6[(2,0)]['real'].get((Fr(-1),Fr(-1)))
print('recorded Res(-a):', rec)
print('-3*Gm(-a) =', sp.simplify(-3*tot))
