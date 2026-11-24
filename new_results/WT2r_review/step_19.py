code = '''
import sys, pickle
sys.path.insert(0,'.')
from fastred import *
from fractions import Fraction as Fr
record = pickle.load(open('record_W2b.pkl','rb'))

# (V1) fixed-point system (*) vs reducer W at each level
def Gm(k,key):
    rec = record.get((k,key))
    return rec['Gm'] if rec and 'Gm' in rec else None
bad=[]; checked=0
for k in range(2,13):
    for key in [kk for kk in record if kk[0]==k]:
        _,(i,j)=key
        if i+j<1: continue
        W = record[key]['W']
        rhs = record[(i+j,(i,j))]['W'] if k==i+j else pf_zero()
        g = Gm(k,(i+1,j))
        if g is not None and i+j+1<=k:
            rhs = pf_add(rhs, pf_scale(pf_mul_DL(g,DL1),Fr(-(i+1))))
        if i>j:
            g = Gm(k,(i,j+1))
            if g is not None:
                rhs = pf_add(rhs, pf_scale(pf_mul_DL(g,DL2),Fr(-(j+1))))
        if j==i-1:
            g = Gm(k,(i,i))
            if g is not None:
                rhs = pf_add(rhs, pf_scale(pf_fold(pf_mul_DL(g,DL1)),Fr(-i)))
        checked+=1
        dff = pf_sub(W,rhs)
        nonzero = (dff[0] and any(not c.is_zero() for c in dff[0].values())) or any(not e.is_zero() for e in dff[1]) or any(not e.is_zero() for e in dff[2]) or not dff[3].is_zero()
        if nonzero: bad.append((k,(i,j)))
print(f"(V1) fixed-point system (*) vs reducer: checked {checked} (k=2..12):", "ALL MATCH" if not bad else bad[:6])

# (V2) pole-order bounds
viol=[]
for (k,(i,j)),rec in record.items():
    d=i+j
    H=rec['H']
    if len(H[1])>2*d-2: viol.append(('W',k,(i,j),len(H[1])))
    if 'Gm' in rec:
        G=rec['Gm']
        if len(G[1])>2*d-3: viol.append(('Gm',k,(i,j),len(G[1])))
print("(V2) pole-order bounds HigherPart(W)<=2d-2, Gm<=2d-3, all 239 records:", "ALL HOLD" if not viol else viol[:8])

# Also: confirm HigherPart has NO order-1 piece by construction and no real poles (pure ipart)
bad2=[]
for (k,key),rec in record.items():
    H=rec['H']
    if H[0] or not H[3].is_zero(): bad2.append((k,key,'H has real/const'))
    if H[1] and not H[1][0].is_zero(): bad2.append((k,key,'H has order-1 at +i'))
print("(V2b) HigherPart is pure ipart with zero order-1 term:", "ALL OK" if not bad2 else bad2[:5])
'''
open('rev_01.py','w').write(code)
import subprocess
print(subprocess.run([sys.executable,'rev_01.py'],capture_output=True,text=True,cwd='.').stdout)
