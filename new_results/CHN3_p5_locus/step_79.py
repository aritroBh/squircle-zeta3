import sys; sys.path.insert(0,'.')
import numpy as np, importlib
import modengine as ME
importlib.reload(ME)
p=2147483647
a=(1,1,1,1,1)
a_arrays=tuple(np.array([a[i]],dtype=np.int64) for i in range(5))
# replicate volreg_mod internals to capture Est, Inn, Gg
import mzv_exact as M
def debug_inn(a_arrays, Wmax, p):
    a1,a2,a3,a4,a5 = a_arrays
    S4=(a1+a2+a3+a4)%p; X24=(a2+a4)%p
    invp=[0]*(Wmax+2)
    for j in range(1,Wmax+2): invp[j]=pow(j,p-2,p)
    def powers(x):
        out={0:np.ones_like(x),1:x}
        for k in range(2,Wmax+1): out[k]=(out[k-1]*x)%p
        return out
    Pa5=powers(a5); Pa4=powers(a4); PX=powers(X24); PS=powers(S4)
    def expH(lc):
        E=[{} for _ in range(Wmax+1)]; E[0]={():np.ones_like(a1)}
        for j in range(1,Wmax+1):
            acc={}
            for l in range(1,j+1):
                scl=(lc[l]*invp[j])%p
                for w,q in ME.st_letter(l,E[j-l],p).items():
                    v=(q*scl)%p
                    if w in acc: acc[w]=(acc[w]+v)%p
                    else: acc[w]=v
            E[j]=acc
        return E
    E5=expH([None]+[((-Pa5[l])%p) for l in range(1,Wmax+1)])
    E24=expH([None]+[((((-1)**(l+1))%p)*PX[l])%p for l in range(1,Wmax+1)])
    ES=expH([None]+[((((-1)**l)%p)*PS[l])%p for l in range(1,Wmax+1)])
    Est=[{} for _ in range(Wmax+1)]; Est[0]={():np.ones_like(a1)}
    for j in range(1,Wmax+1):
        acc={}
        for j5 in range(j+1):
            for j24 in range(j-j5+1):
                jS=j-j5-j24
                d=ME.st_dicts(E5[j5],E24[j24],p); d=ME.st_dicts(d,ES[jS],p)
                for w,x in d.items():
                    if w in acc: acc[w]=(acc[w]+x)%p
                    else: acc[w]=x
        Est[j]=acc
    return Est, E5, E24, ES
Est,E5,E24,ES = debug_inn(a_arrays,5,p)
print("Est[2]:", {w:int(v[0]) for w,v in Est[2].items()}, " expected {(1,1):9,(2,):10}")
print("E5[1]:", {w:int(v[0]) for w,v in E5[1].items()}, " E24[1]:", {w:int(v[0]) for w,v in E24[1].items()}, " ES[1]:", {w:int(v[0]) for w,v in ES[1].items()})
print("E5[2]:", {w:int(v[0]) for w,v in E5[2].items()})
print("E24[2]:", {w:int(v[0]) for w,v in E24[2].items()})
print("ES[2]:", {w:int(v[0]) for w,v in ES[2].items()})
