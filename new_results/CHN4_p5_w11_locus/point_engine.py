"""Scalar point engine for P5: reduced DS coords at a rational point a=(a1..a5), W<=Wmax.
Mirrors engine5.build_engine/volreg_words but with Fraction scalars instead of poly dicts."""
from fractions import Fraction as F
import mzv_exact as M

_TABLES={}
def table(W):
    if W not in _TABLES:
        _TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return _TABLES[W]

def st_letter(i, d):
    nd={}
    for w,c in d.items():
        for w2,c2 in M.stuffle((i,),w).items():
            nd[w2]=nd.get(w2,F(0))+c*c2
    return {w:c for w,c in nd.items() if c}
def st_dicts(d1,d2):
    nd={}
    for w1,c1 in d1.items():
        for w2,c2 in d2.items():
            for w,c in M.stuffle(w1,w2).items():
                nd[w]=nd.get(w,F(0))+c1*c2*c
    return {w:c for w,c in nd.items() if c}
def prepend(c, d):
    return {(c,)+w: p for w,p in d.items()}

def volreg_words_point(a, Wmax):
    a1,a2,a3,a4,a5 = [F(x) for x in a]
    S4 = a1+a2+a3+a4
    X24 = a2+a4
    def expH(lc, Wmax):  # lc[l] scalar
        E=[{} for _ in range(Wmax+1)]; E[0]={():F(1)}
        for j in range(1,Wmax+1):
            acc={}
            for l in range(1,j+1):
                d=st_letter(l, E[j-l])
                for w,p in d.items():
                    acc[w]=acc.get(w,F(0))+p*lc[l]/j
            E[j]={w:p for w,p in acc.items() if p}
        return E
    E5 = expH([None]+[-a5**l for l in range(1,Wmax+1)], Wmax)
    E24= expH([None]+[((-1)**(l+1))*X24**l for l in range(1,Wmax+1)], Wmax)
    ES = expH([None]+[((-1)**l)*S4**l for l in range(1,Wmax+1)], Wmax)
    # p_n (n>=2), q_m
    p={}
    for n in range(2,Wmax+1):
        p[n]=((-1)**n)*a4*(X24**(n-1)-a4**(n-1))
    q={0:F(1)}
    for m in range(1,Wmax+1):
        q[m]=((-1)**m)*(S4**m - X24*S4**(m-1))
    r={}
    for cm1 in range(2,Wmax):
        acc=F(0)
        for n in range(2,cm1+1):
            acc+=p[n]*q[cm1-n]
        r[cm1]=-a5*acc
    Est=[{} for _ in range(Wmax+1)]; Est[0]={():F(1)}
    for j in range(1,Wmax+1):
        acc={}
        for j5 in range(j+1):
            for j24 in range(j-j5+1):
                jS=j-j5-j24
                d=st_dicts(E5[j5],E24[j24]); d=st_dicts(d,ES[jS])
                for w,pp in d.items(): acc[w]=acc.get(w,F(0))+pp
        Est[j]={w:pp for w,pp in acc.items() if pp}
    Inn=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for c in range(3,W+1):
            if (c-1) not in r: continue
            d=prepend(c, Est[W-c])
            for w,pp in d.items(): acc[w]=acc.get(w,F(0))+r[c-1]*pp
        Inn[W]={w:pp for w,pp in acc.items() if pp}
    A13=a1+a3
    Gg=[{} for _ in range(Wmax+1)]; Gg[0]={():F(1)}
    for j in range(1,Wmax+1):
        acc={}
        for m in range(2,j+1):
            lm=((-1)**m)*(A13**m + X24**m - S4**m)
            d=st_letter(m, Gg[j-m])
            for w,pp in d.items(): acc[w]=acc.get(w,F(0))+pp*lm/j
        Gg[j]={w:pp for w,pp in acc.items() if pp}
    Areg=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for i in range(0,W+1):
            j=W-i
            if not Gg[i] or not Inn[j]: continue
            d=st_dicts(Gg[i],Inn[j])
            for w,pp in d.items(): acc[w]=acc.get(w,F(0))+pp
        Areg[W]={w:pp for w,pp in acc.items() if pp}
    # volreg: Areg + reversed
    out=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        out[W]=dict(Areg[W])
    return out

def reduced_at(a, Wmax, forms_words=None):
    """full reduced coords at point a: regular (A+A_rev) + m0."""
    vw = volreg_words_point(a, Wmax)
    vwr = volreg_words_point(a[::-1], Wmax)   # reversed point gives A_rev
    import kimilib
    m0 = kimilib.m0_coords(a, Wmax)
    out={}
    for W in range(2,Wmax+1):
        T=table(W); d=T.nullity
        acc=[F(0)]*d
        for w,pp in vw[W].items():
            rw=T.reduce_word(w)
            for b in range(d): acc[b]+=pp*rw[b]
        for w,pp in vwr[W].items():
            rw=T.reduce_word(w)
            for b in range(d): acc[b]+=pp*rw[b]
        for b in range(d): acc[b]+=m0[W][b]
        out[W]=acc
    return out
