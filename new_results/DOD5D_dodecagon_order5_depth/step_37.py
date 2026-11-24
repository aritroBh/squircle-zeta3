import sys
sys.path.insert(0, ".")
import importlib, mpmath as mp
mp.mp.dps = 50
import dod5d_exact as DX
importlib.reload(DX)
r3 = mp.sqrt(3)
KAPm = [mp.mpf(str(k[0]))+mp.mpf(str(k[1]))*r3 for k in DX.KAP]
C0m  = [mp.mpf(str(c[0]))+mp.mpf(str(c[1]))*r3 for c in DX.C0S]
print("C0S[2] should be (2-r3)/2:", C0m[2], (2-r3)/2)
print("KAP:", KAPm)
print("LETTERS:", [mp.mpf(str(l[0]))+mp.mpf(str(l[1]))*r3 for l in DX.LETTERS])
def eval_d(d, t):
    LT=mp.log(t); LU=mp.log(1-t)
    W=[mp.log(1+KAPm[i]*t) for i in range(4)]
    CA=mp.log((2-r3)/2); CC0=mp.log(2-r3); CC1=mp.log((2*r3-3)/2)
    logs=[LT,LU,W[0],W[1],W[2],W[3],CA,CC0,CC1]
    tot=mp.mpf(0)
    for (n,le),c in d.items():
        v=t**n
        for i,e in enumerate(le):
            if e: v*=logs[i]**e
        tot+=(mp.mpf(str(c[0]))+mp.mpf(str(c[1]))*r3)*v
    return tot
t = mp.mpf('0.37')
ap=[C0m[i]*KAPm[i] for i in range(4)]
CCv = {7: mp.log(2-r3), 8: mp.log((2*r3-3)/2), 6: mp.log((2-r3)/2)}
la=[CCv[DX.CCslot[i]]+mp.log(1+KAPm[i]*t) for i in range(4)]
d = DX.taugam('gam',1,1)
print("dict gam11:", mp.nstr(eval_d(d,t),25))
print("ref  gam11:", mp.nstr(sum(ap[i]*la[i] for i in range(4)),25))
# inspect raw dict
for k,v in sorted(d.items()):
    print(k, v)
