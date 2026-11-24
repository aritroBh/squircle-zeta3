modeng = '''"""Vectorized modular point engine for P5: c_W(a) mod p at K points simultaneously.
Scalars: numpy int64 arrays (length K), arithmetic mod p. Exact via rational reconstruction."""
import numpy as np
from fractions import Fraction as F
import mzv_exact as M

_TABLES={}
def table(W):
    if W not in _TABLES:
        _TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return _TABLES[W]

_RWCACHE={}
def reduce_word_mod(W, w, p):
    """reduce_word of w as mod-p integers (length d)."""
    key=(W,w)
    c=_RWCACHE.get(key)
    if c is None:
        T=table(W)
        rw=T.reduce_word(w)  # Fractions
        _RWCACHE[key]=rw
        c=rw
    return [ (x.numerator%p)*pow(x.denominator%p, p-2, p)%p if x else 0 for x in c]

def st_letter(i, d, p):
    nd={}
    for w,c in d.items():
        for w2,c2 in M.stuffle((i,),w).items():
            if w2 in nd: nd[w2]=(nd[w2]+c*c2)%p
            else: nd[w2]=(c*c2)%p
    return nd
def st_dicts(d1,d2,p):
    nd={}
    for w1,c1 in d1.items():
        for w2,c2 in d2.items():
            for w,c in M.stuffle(w1,w2).items():
                if w in nd: nd[w]=(nd[w]+c1*c2)%p
                else: nd[w]=(c1*c2)%p
    return nd

def _ones(a): return np.ones_like(a)

def volreg_mod(a_arrays, Wmax, p):
    a1,a2,a3,a4,a5 = a_arrays
    S4=(a1+a2+a3+a4)%p; X24=(a2+a4)%p
    invp=[0]*(Wmax+2)
    for j in range(1,Wmax+2): invp[j]=pow(j,p-2,p)
    def powers(x):
        out={0:_ones(x),1:x}
        for k in range(2,Wmax+1): out[k]=(out[k-1]*x)%p
        return out
    Pa5=powers(a5); Pa4=powers(a4); PX=powers(X24); PS=powers(S4); PA13=powers((a1+a3)%p)
    def expH(lc):
        E=[{} for _ in range(Wmax+1)]; E[0]={():_ones(a1)}
        for j in range(1,Wmax+1):
            acc={}
            for l in range(1,j+1):
                scl=(lc[l]*invp[j])%p
                for w,q in st_letter(l,E[j-l],p).items():
                    v=(q*scl)%p
                    if w in acc: acc[w]=(acc[w]+v)%p
                    else: acc[w]=v
            E[j]=acc
        return E
    E5=expH([None]+[((-Pa5[l])%p) for l in range(1,Wmax+1)])
    E24=expH([None]+[((((-1)**(l+1))%p)*PX[l])%p for l in range(1,Wmax+1)])
    ES=expH([None]+[((((-1)**l)%p)*PS[l])%p for l in range(1,Wmax+1)])
    pp={}
    for n in range(2,Wmax+1):
        pp[n]=((((-1)**n)%p)*a4)%p * ((PX[n-1]-Pa4[n-1])%p) % p
    qq={0:_ones(a1)}
    for m in range(1,Wmax+1):
        qq[m]=((((-1)**m)%p)*((PS[m]-(X24*PS[m-1])%p)%p))%p
    rr={}
    for cm1 in range(2,Wmax):
        acc=np.zeros_like(a1)
        for n in range(2,cm1+1): acc=(acc+pp[n]*qq[cm1-n])%p
        rr[cm1]=((p-a5)%p)*acc%p
    Est=[{} for _ in range(Wmax+1)]; Est[0]={():_ones(a1)}
    for j in range(1,Wmax+1):
        acc={}
        for j5 in range(j+1):
            for j24 in range(j-j5+1):
                jS=j-j5-j24
                d=st_dicts(E5[j5],E24[j24],p); d=st_dicts(d,ES[jS],p)
                for w,x in d.items():
                    if w in acc: acc[w]=(acc[w]+x)%p
                    else: acc[w]=x
        Est[j]=acc
    Inn=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for c in range(3,W+1):
            if (c-1) not in rr: continue
            for w,x in (( (c,)+w, x) for w,x in Est[W-c].items()):
                v=(x*rr[c-1])%p
                if w in acc: acc[w]=(acc[w]+v)%p
                else: acc[w]=v
        Inn[W]=acc
    Gg=[{} for _ in range(Wmax+1)]; Gg[0]={():_ones(a1)}
    for j in range(1,Wmax+1):
        acc={}
        for m in range(2,j+1):
            lm=(((-1)**m)%p)*((PA13[m]+PX[m]-PS[m])%p)%p
            scl=(lm*invp[j])%p
            for w,x in st_letter(m,Gg[j-m],p).items():
                v=(x*scl)%p
                if w in acc: acc[w]=(acc[w]+v)%p
                else: acc[w]=v
        Gg[j]=acc
    Areg=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for i in range(0,W+1):
            j=W-i
            if not Gg[i] or not Inn[j]: continue
            for w,x in st_dicts(Gg[i],Inn[j],p).items():
                if w in acc: acc[w]=(acc[w]+x)%p
                else: acc[w]=x
        Areg[W]=acc
    return Areg

def gamprod_series_mod(coeffs, Wmax, p):
    """coeffs: list of (ell_array, e) with e in {+1,-1}. q[m]=(-1)^m/m sum e*ell^m (array)."""
    invp={j:pow(j,p-2,p) for j in range(1,Wmax+1)}
    q={}
    for m in range(2,Wmax+1):
        s=np.zeros_like(coeffs[0][0])
        for ell,e in coeffs:
            pw=ell.copy()
            for _ in range(m-1): pw=(pw*ell)%p
            s=(s+ (e%p)*pw)%p
        q[m]=((((-1)**m)%p)*invp[m]%p)*s%p
    E=[{} for _ in range(Wmax+1)]; E[0]={():_ones(coeffs[0][0])}
    for W in range(1,Wmax+1):
        acc={}
        for i in range(2,W+1):
            for z,c in q[i].items():
                for z2,c2 in E[W-i].items():
                    z3=tuple(sorted(z+z2))
                    v=(c*c2)%p*(i*invp[W]%p)%p
                    if z3 in acc: acc[z3]=(acc[z3]+v)%p
                    else: acc[z3]=v
        E[W]=acc
    return E

def m0_coords_mod(a_arrays, Wmax, p):
    a1,a2,a3,a4,a5=a_arrays
    S4=(a1+a2+a3+a4)%p; S4r=(a5+a4+a3+a2)%p
    s1=gamprod_series_mod([( (a2+a4)%p ,1),(S4,-1),((a1+a3)%p,1)],Wmax,p)
    s2=gamprod_series_mod([((a2+a4)%p,1),(S4r,-1),((a5+a3)%p,1)],Wmax,p)
    inv24=np.array([pow(int(v),p-2,p) for v in (a2+a4)%p],dtype=np.int64)
    f1=(a2*inv24)%p; f2=(a4*inv24)%p
    out={}
    for W in range(2,Wmax+1):
        d=table(W).nullity
        acc=[np.zeros_like(a1) for _ in range(d)]
        for s_,f_ in ((s1,f1),(s2,f2)):
            for zmon,c in s_[W].items():
                cc=(c*f_)%p
                wc=M.wordcombo(tuple((k,) for k in zmon))
                for w,cw in wc.items():
                    rw=reduce_word_mod(W,w,p)
                    cwmod=(cw.numerator%p)*pow(cw.denominator%p,p-2,p)%p
                    for b in range(d):
                        if rw[b]:
                            acc[b]=(acc[b]+cc*cwmod*rw[b])%p
        out[W]=acc
    return out

def full_coords_mod(a_arrays, Wmax, p):
    """Returns {W: [basis coord arrays mod p]} = regular (A+A_rev) + m0."""
    arev=tuple(a_arrays[::-1])
    vw=volreg_mod(a_arrays,Wmax,p); vwr=volreg_mod(arev,Wmax,p)
    m0=m0_coords_mod(a_arrays,Wmax,p)
    out={}
    for W in range(2,Wmax+1):
        d=table(W).nullity
        acc=[np.zeros_like(a_arrays[0]) for _ in range(d)]
        if W>=3:
            for src in (vw[W], vwr[W]):
                for w,x in src.items():
                    rw=reduce_word_mod(W,w,p)
                    for b in range(d):
                        if rw[b]: acc[b]=(acc[b]+x*rw[b])%p
        for b in range(d):
            acc[b]=(acc[b]+m0[W][b])%p
        out[W]=acc
    return out
'''
open('modengine.py','w').write(modeng)
print("modengine.py rewritten cleanly")
