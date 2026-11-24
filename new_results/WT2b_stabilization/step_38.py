import sys, pickle
sys.path.insert(0,'.')
from fastred import *
record = pickle.load(open('record_W2b.pkl','rb'))

# (V1) verify the triangular fixed-point system (*) reproduces the reducer's W at each level:
# W_ij = delta_{k,d} R_ij - (i+1) Gm_{i+1,j} DL1 - [i>j](j+1) Gm_{i,j+1} DL2 - [j=i-1] i F(Gm_ii DL1)
# where R_ij = W_ij at k = i+j (initial level, feeds vanish there)
def Gm(k,key):
    rec = record.get((k,key))
    return rec['Gm'] if rec and 'Gm' in rec else None
bad=[]
checked=0
for k in range(2,13):
    dmax=k
    for key in [kk for kk in record if kk[0]==k]:
        _,(i,j)=key
        if i+j<1: continue
        W = record[key]['W']
        rhs = record[(i+j,(i,j))]['W'] if k==i+j else pf_zero()
        # wait: W_ij at k=i+j is R_ij only if feeds vanish: at level d, feed keys have degree d+1>d: absent. OK
        g = Gm(k,(i+1,j))
        if g is not None and i+j+1<=k:
            t = pf_scale(pf_mul_DL(g,DL1),Fr(-(i+1)))
            rhs = pf_add(rhs,t)
        if i>j:
            g = Gm(k,(i,j+1))
            if g is not None:
                rhs = pf_add(rhs, pf_scale(pf_mul_DL(g,DL2),Fr(-(j+1))))
        if j==i-1:
            g = Gm(k,(i,i))
            if g is not None:
                rhs = pf_add(rhs, pf_scale(pf_fold(pf_mul_DL(g,DL1)),Fr(-i)))
        checked+=1
        # compare
        dff = pf_sub(W,rhs)
        nonzero = (dff[0] and any(not c.is_zero() for c in dff[0].values())) or any(not e.is_zero() for e in dff[1]) or any(not e.is_zero() for e in dff[2]) or not dff[3].is_zero()
        if nonzero: bad.append((k,(i,j)))
print(f"(V1) fixed-point system (*) vs reducer: checked {checked} (k=2..12):", "ALL MATCH" if not bad else bad[:6])

# (V2) pole-order bounds: HigherPart(W_ij) <= 2d-2, Gm_ij <= 2d-3
viol=[]
for (k,(i,j)),rec in record.items():
    d=i+j
    H=rec['H']
    if len(H[1])>2*d-2: viol.append(('W',k,(i,j),len(H[1])))
    if 'Gm' in rec:
        G=rec['Gm']
        if len(G[1])>2*d-3: viol.append(('Gm',k,(i,j),len(G[1])))
print("(V2) pole-order bounds on all data:", "ALL HOLD" if not viol else viol[:8])

# (V3a) top-order recursion: h_ij(k) for k>d satisfies the recursion with tau1_0,tau2_0
tau1 = QSI(-1,0,-1,0)/QSI(2); tau2 = QSI(-1,0,1,0)/QSI(2)
def h(k,key):
    i,j=key; d=i+j
    rec=record.get((k,key))
    if rec is None: return None
    H=rec['H'][1]; top=2*d-2
    return H[top-1] if len(H)>=top else ZERO
bad=[]
for k in range(3,13):
    for i in range(1,k+1):
        for j in range(0,i+1):
            d=i+j
            if d<2 or d>k: continue
            if k==d: continue
            hv = h(k,(i,j))
            if hv is None: continue
            # recursion: h = [(i+1)tau1 h_{i+1,j} + [i>j](j+1)tau2 h_{i,j+1} + [j=i-1] i tau2 (h_ii - conj h_ii)]/(2d-1)
            s = ZERO
            hA = h(k,(i+1,j))
            if hA is not None: s = s + hA*QSI(i+1)*tau1
            if i>j:
                hB = h(k,(i,j+1))
                if hB is not None: s = s + hB*QSI(j+1)*tau2
            if j==i-1:
                hd = h(k,(i,i))
                if hd is not None: s = s + (hd - hd.conj_i())*QSI(i)*tau2
            s = s/QSI(2*d-1)
            if hv != s: bad.append((k,(i,j),hv.sympy(),s.sympy()))
print("(V3a) top-order recursion vs data (k=3..12):", "ALL MATCH" if not bad else bad[:6])
