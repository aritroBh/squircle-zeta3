import os, signal
try: os.kill(int(open('eng12b.pid').read()), signal.SIGKILL); print("killed build")
except Exception as e: print(e)
modeng = '''"""Vectorized modular point engine for P5: c_W(a) mod p at K points simultaneously.
Scalars are numpy int64 arrays (length K), arithmetic mod p. Exact via rational reconstruction.
Mirrors engine5.build_engine/volreg_words + kimilib.m0_coords."""
import numpy as np
from fractions import Fraction as F
import mzv_exact as M

class ModCtx:
    def __init__(self, p):
        self.p=int(p)
    def inv_tab(self, n):
        p=self.p
        inv=np.empty(n+1, dtype=np.int64)
        for j in range(1,n+1): inv[j]=pow(j, p-2, p)
        return inv

_TABLES={}
def table(W):
    if W not in _TABLES:
        _TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return _TABLES[W]

def st_letter(i, d):
    nd={}
    for w,c in d.items():
        for w2,c2 in M.stuffle((i,),w).items():
            v=c*c2
            k=w2
            if k in nd: nd[k]=(nd[k]+v)%MODP
            else: nd[k]=v
    return {w:c for w,c in nd.items()}
def st_dicts(d1,d2):
    nd={}
    for w1,c1 in d1.items():
        for w2,c2 in d2.items():
            for w,c in M.stuffle(w1,w2).items():
                v=(c1*c2)%MODP
                if w in nd: nd[w]=(nd[w]+v)%MODP
                else: nd[w]=v
    return nd
def prepend(c, d):
    return {(c,)+w: p for w,p in d.items()}

MODP=None
def volreg_words_point_mod(a_arrays, Wmax, p):
    """a_arrays: tuple of 5 int64 arrays (mod p). Returns word->array dicts per weight."""
    global MODP; MODP=p
    a1,a2,a3,a4,a5 = a_arrays
    S4 = (a1+a2+a3+a4)%p
    X24 = (a2+a4)%p
    powtab = ModCtx(p).inv_tab(Wmax+1)   # modular inverses of 1..Wmax+1
    def expH(lc):  # lc[l]: array
        E=[{} for _ in range(Wmax+1)]; E[0]={():None}  # placeholder
        E[0]={(): np.ones_like(a1)}
        for j in range(1,Wmax+1):
            acc={}
            for l in range(1,j+1):
                d=st_letter(l, E[j-l])
                scl = (lc[l]*powtab[j])%p
                for w,q in d.items():
                    v=(q*scl)%p
                    if w in acc: acc[w]=(acc[w]+v)%p
                    else: acc[w]=v
            E[j]=acc
        return E
    powa5={}; powX24={}; powS4={}
    pw=a1.copy(); 
    # powers
    def powers(x):
        out={1:x}
        for k in range(2,Wmax+1): out[k]=(out[k-1]*x)%p
        return out
    Pa5=powers(a5); PX24=powers(X24); PS4=powers(S4)
    E5 = expH([None]+[((-1)*Pa5[l])%p for l in range(1,Wmax+1)])
    E24= expH([None]+[(((-1)**(l+1))*PX24[l])%p for l in range(1,Wmax+1)])
    ES = expH([None]+[(((-1)**l)*PS4[l])%p for l in range(1,Wmax+1)])
    p_={}
    for n in range(2,Wmax+1):
        p_[n]=(((-1)**n)*a4%p)*((PX24[n-1]-Pa5 if False else 0))  # placeholder fix below
    # p_n = (-1)^n a4 (X24^{n-1} - a4^{n-1})
    Pa4=powers(a4)
    p_={}
    for n in range(2,Wmax+1):
        p_[n]=((((-1)**n)%p)*a4*( (PX24[n-1]-Pa4[n-1])%p ))%p
    q={0:np.ones_like(a1)}
    for m in range(1,Wmax+1):
        q[m]=(((-1)**m)%p)*((PS4[m]-(X24*PS4[m-1])%p)%p)%p
    r={}
    for cm1 in range(2,Wmax):
        acc=np.zeros_like(a1)
        for n in range(2,cm1+1):
            acc=(acc+p_[n]*q[cm1-n])%p
        r[cm1]=((-a5)%p)*acc%p
    Est=[{} for _ in range(Wmax+1)]; Est[0]={():np.ones_like(a1)}
    for j in range(1,Wmax+1):
        acc={}
        for j5 in range(j+1):
            for j24 in range(j-j5+1):
                jS=j-j5-j24
                d=st_dicts(E5[j5],E24[j24]); d=st_dicts(d,ES[jS])
                for w,pp in d.items():
                    if w in acc: acc[w]=(acc[w]+pp)%p
                    else: acc[w]=pp
        Est[j]=acc
    Inn=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for c in range(3,W+1):
            if (c-1) not in r: continue
            d=prepend(c, Est[W-c])
            for w,pp in d.items():
                v=(pp*r[c-1])%p
                if w in acc: acc[w]=(acc[w]+v)%p
                else: acc[w]=v
        Inn[W]=acc
    A13=(a1+a3)%p
    PA13=powers(A13)
    Gg=[{} for _ in range(Wmax+1)]; Gg[0]={():np.ones_like(a1)}
    for j in range(1,Wmax+1):
        acc={}
        for m in range(2,j+1):
            lm=(((-1)**m)%p)*((PA13[m]+PX24[m]-PS4[m])%p)%p
            d=st_letter(m, Gg[j-m])
            scl=(lm*powtab[j])%p
            for w,pp in d.items():
                v=(pp*scl)%p
                if w in acc: acc[w]=(acc[w]+v)%p
                else: acc[w]=v
        Gg[j]=acc
    Areg=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for i in range(0,W+1):
            j=W-i
            if not Gg[i] or not Inn[j]: continue
            d=st_dicts(Gg[i],Inn[j])
            for w,pp in d.items():
                if w in acc: acc[w]=(acc[w]+pp)%p
                else: acc[w]=pp
        Areg[W]=acc
    return {W:Areg[W] for W in range(3,Wmax+1)}
'''
open('modengine.py','w').write(modeng)
print("modengine.py written (part 1)")
