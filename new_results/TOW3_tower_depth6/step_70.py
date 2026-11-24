import sys, os, pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
tn = f'{base}/TOWN1_tower_large_n'
wts = f'{base}/WTS_weight_k_structure/theory'
sys.path.insert(0, tn); sys.path.insert(0, wts)
os.chdir(tn)
import kernel_from_tower as KFT
from fastred import QSI
import sympy as sp
s2 = sp.sqrt(2); ii = sp.I
av = 1+s2
def qsi_sp(q): return sp.Rational(q.p0) + sp.Rational(q.p1)*s2 + sp.Rational(q.p2)*ii + sp.Rational(q.p3)*s2*ii
T6 = pickle.load(open(f'{wts}/tower_6.pkl','rb'))
rec = T6.get((3,0), {})
h = {r: QSI(*rec[r]) for r in rec if r >= 2}
gm = KFT.gm_from_tower(h)
tot = qsi_sp(gm[3])
for idx, e in enumerate(gm[1]):
    if idx == 0: continue
    tot += qsi_sp(e) * (-av-ii)**(-idx)
for idx, e in enumerate(gm[2]):
    if idx == 0: continue
    tot += qsi_sp(e) * (-av+ii)**(-idx)
tot = sp.expand(tot)
K6 = pickle.load(open(f'{wts}/kern_6.pkl','rb'))
rec = K6[(2,0)]['real'].get((Fr(-1),Fr(-1)))
print('Gm(-a):', sp.simplify(tot))
print('-3*Gm(-a):', sp.simplify(-3*tot), '  recorded:', rec)
