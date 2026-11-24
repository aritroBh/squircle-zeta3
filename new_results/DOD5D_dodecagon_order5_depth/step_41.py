import pickle, sys
sys.path.insert(0, "<project>/research/08_theory/E5c_nonlocal_z5")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 120
import symeng, dod5d_exact as DX
E5X = pickle.load(open("E5_dod_expanded.pkl","rb"))
r3 = mp.sqrt(3)
M = []
for j in range(6):
    c, s = mp.cos(j*mp.pi/6), mp.sin(j*mp.pi/6)
    M.append(((1 if j<=3 else -1)*c, (1 if j<=3 else -1)*s))
y0 = (2-r3, mp.mpf(0)); y1 = ((2*r3-3)/2, (2-r3)/2)
P = symeng.edge_params(M, y0, y1)
u = mp.sqrt(2-r3)
KAPm = [mp.mpf(str(k[0]))+mp.mpf(str(k[1]))*r3 for k in DX.KAP]
logCA = mp.log((2-r3)/2); logCC0 = mp.log(2-r3); logCC1 = mp.log((2*r3-3)/2)
def E5_exact_eval(t):
    LT = mp.log(t); LU = mp.log(1-t)
    W = [mp.log(1+KAPm[i]*t) for i in range(4)]
    tot = mp.mpf(0)
    for (n, le), (p,q) in E5X.items():
        v = t**n * LT**le[0] * LU**le[1]
        for i in range(4):
            if le[2+i]: v *= W[i]**le[2+i]
        if le[6]: v *= logCA**le[6]
        if le[7]: v *= logCC0**le[7]
        if le[8]: v *= logCC1**le[8]
        tot += (mp.mpf(str(p))+mp.mpf(str(q))*r3)*v
    return u*tot
for tt in ['0.37','0.71']:
    t = mp.mpf(tt)
    a = E5_exact_eval(t); b = symeng.E5_sym(t, M, P, conv='B')
    print(f"t={tt}: a={mp.nstr(a,30)}")
    print(f"        b={mp.nstr(b,30)}  diff={mp.nstr(abs(a-b),3)}")
