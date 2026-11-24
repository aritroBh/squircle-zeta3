import sys, os, pickle
from fractions import Fraction as Fr
base = '<project>/research/08_theory'
tn = f'{base}/TOWN1_tower_large_n'; wts = f'{base}/WTS_weight_k_structure/theory'
sys.path.insert(0, tn); sys.path.insert(0, wts); os.chdir(tn)
import kernel_from_tower as KFT
from fastred import QSI, ZERO
T6 = pickle.load(open(f'{wts}/tower_6.pkl','rb'))
def hd(key):
    rec = T6.get(key, {})
    return {r: QSI(*rec[r]) for r in rec if r >= 2}
K20, K11, H20, H11 = KFT.kernel_20_21_level(hd((3,0)), hd((2,1)))
K6 = pickle.load(open(f'{wts}/kern_6.pkl','rb'))
rec = K6[(2,0)]['real']
print('recorded -a:', rec.get((Fr(-1),Fr(-1))))
print('computed K20 real at -a:', K20[0].get(-1-__import__('fastred').SQ2 if hasattr(__import__('fastred'),'SQ2') else None, 'n/a'))
# print all computed real poles of K20
from fastred import AQ, BQ
for r, c in K20[0].items():
    print('computed pole', (r.p0, r.p1), ':', (c.p0, c.p1, c.p2, c.p3))
print('AQ,BQ =', (AQ.p0,AQ.p1),(BQ.p0,BQ.p1))
