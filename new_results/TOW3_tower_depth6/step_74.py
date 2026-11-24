import sys, os, pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
tn = f'{base}/TOWN1_tower_large_n'; wts = f'{base}/WTS_weight_k_structure/theory'
sys.path.insert(0, tn); sys.path.insert(0, wts); os.chdir(tn)
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
print('const:', (gm[3].p0, gm[3].p1, gm[3].p2, gm[3].p3))
print('ipos:', [(e.p0,e.p1,e.p2,e.p3) for e in gm[1]])
print('ineg:', [(e.p0,e.p1,e.p2,e.p3) for e in gm[2]])
# residue at -a per pf_mul_DL: const*1 + sum_j e_j/(-a-i)^j + f_j/(-a+i)^j
tot = qsi_sp(gm[3])
for j, e in enumerate(gm[1], start=1):
    tot += qsi_sp(e)/(-av-ii)**j
for j, e in enumerate(gm[2], start=1):
    tot += qsi_sp(e)/(-av+ii)**j
tot = sp.expand(tot)
print('Res(-a) of Gm DL1:', sp.simplify(tot))
print('times -3:', sp.simplify(-3*tot))
