
with open('progress.log','a') as f:
    f.write("2026-10-29 PartC formulas: t-integral == direct x-integral to 51 digits (pt1) and >72 digits (pt2, dps80); A+A_rev single series 22 digits; double series exact at anchor u=(1..1) (=2/15) and ~12 digits at both pts (Richardson tail). Two self-corrected index bugs (S4->S3 in K; K off-by-one). numcheck_pt1/pt2/dbl.pkl\n")

# Now write engine5.py
code = '''"""Exact symbolic engine for the 5-chain P5: c_W(a1..a5) = [u^W] Vol/32 as
DS-basis coordinate polys (dicts (e1..e5)->Fraction), for W <= Wmax.
Vol/32 = A(u1,u2;u4,u5) + A(u5,u4;u2,u1),
A = G~ sum_{m>=1} alpha_m  +  (m=0 term, handled separately: Gamma-product, no irreducibles)
G~ = Gamma(1+u1+u3)Gamma(1+u2+u4)/Gamma(1+u1+u2+u3+u4)
alpha_m = (-u5)/m P5_{m-1} [u4u2/((u4+m)(u2+u4+m))] J_m P24_{m-1} PS_{m-1}
P5 = prod(1-u5/j), P24 = prod(1+(u2+u4)/j), PS = prod(1+S4u/j)^-1, J_m=(1+(u2+u4)/m)/(1+S4u/m)"""
from fractions import Fraction as F
import mzv_exact as M

def padd(p,q,s=F(1)):
    out=dict(p)
    for e,c in q.items():
        out[e]=out.get(e,F(0))+s*c
        if out[e]==0: del out[e]
    return out
def pscale(p,s):
    return {} if s==0 else {e:c*s for e,c in p.items()}
def pmul(p,q):
    out={}
    for e1,c1 in p.items():
        for e2,c2 in q.items():
            k=tuple(a+b for a,b in zip(e1,e2))
            out[k]=out.get(k,F(0))+c1*c2
    return out
def var(i,k):
    e=[0]*5; e[i]=k; return {tuple(e):F(1)}
ZERO=(0,0,0,0,0)
def linpow5(terms, m):  # (sum c_i a_i)^m
    from math import factorial
    out={}
    def rec(t, rem, exps, coef):
        if t==len(terms)-1:
            i,c=terms[t]
            e=list(exps); e[i]+=rem
            out[tuple(e)]=out.get(tuple(e),F(0))+coef*(c**rem)
            return
        i,c=terms[t]
        for k in range(rem+1):
            e=list(exps); e[i]+=k
            rec(t+1, rem-k, tuple(e), coef*(c**k)*F(factorial(rem),factorial(k)*factorial(rem-k)))
    rec(0,m,(0,0,0,0,0),F(1))
    return out

def st_letter(i, d):
    nd={}
    for w,c in d.items():
        for w2,c2 in M.stuffle((i,),w).items():
            q=nd.get(w2)
            nd[w2]=padd(q,c,c2) if q is not None else pscale(c,c2)
    return {w:c for w,c in nd.items() if c}
def st_dicts(d1,d2):
    nd={}
    for w1,c1 in d1.items():
        for w2,c2 in d2.items():
            for w,c in M.stuffle(w1,w2).items():
                q=nd.get(w)
                pq=pmul(c1,c2)
                nd[w]=padd(q,pq,c) if q is not None else pscale(pq,c)
    return {w:c for w,c in nd.items() if c}
def prepend(c, d):
    return {(c,)+w: p for w,p in d.items()}

def build_engine(Wmax=10):
    S4 = [(0,1),(1,1),(2,1),(3,1)]   # a1+a2+a3+a4
    X24 = [(1,1),(3,1)]              # a2+a4
    # exp-H series
    def expH(lc, Wmax):
        # lc[l] = l*c_l poly; E[j] = (1/j) sum_l st_letter(l, E[j-l]) * lc[l]
        E=[{} for _ in range(Wmax+1)]; E[0]={(): {ZERO:F(1)}}
        for j in range(1,Wmax+1):
            acc={}
            for l in range(1,j+1):
                d=st_letter(l, E[j-l])
                for w,p in d.items():
                    q=acc.get(w)
                    pq=pscale(pmul(p, lc[l]), F(1,j))
                    acc[w]=padd(q,pq) if q is not None else pq
            E[j]={w:p for w,p in acc.items() if p}
        return E
    E5 = expH([None]+[pscale(var(4,l), F(-1)) for l in range(1,Wmax+1)], Wmax)
    E24= expH([None]+[pscale(linpow5(X24,l), F((-1)**(l+1))) for l in range(1,Wmax+1)], Wmax)
    ES = expH([None]+[pscale(linpow5(S4,l), F((-1)**l)) for l in range(1,Wmax+1)], Wmax)
    # p_n (n>=2), q_m (m>=0)
    p={}
    for n in range(2,Wmax+1):
        pn=padd(linpow5(X24,n-1), var(1,n-1), F(-1))
        p[n]=pscale(pmul(var(3,1), pn), F((-1)**(n+1)))
    q={}
    for m in range(0,Wmax+1):
        acc={}
        for i in range(m+1):
            k=m-i
            term=pmul(linpow5(X24,i), linpow5(S4,k))
            acc=padd(acc, term, F((-1)**k))
        q[m]=acc
    # r_{c-1} = -a5 sum_{n+m=c-1, n>=2, m>=0} p_n q_m
    r={}
    for cm1 in range(2,Wmax):   # c-1 >= 2 (n>=2)
        acc={}
        for n in range(2, cm1+1):
            m=cm1-n
            acc=padd(acc, pmul(p[n], q[m]))
        r[cm1]=pscale(pmul(var(4,1), acc), F(-1))
    # Est[j] = sum_{j5+j24+jS=j} stuffle3
    Est=[{} for _ in range(Wmax+1)]; Est[0]={(): {ZERO:F(1)}}
    for j in range(1,Wmax+1):
        acc={}
        for j5 in range(j+1):
            for j24 in range(j-j5+1):
                jS=j-j5-j24
                d=st_dicts(E5[j5], E24[j24])
                d=st_dicts(d, ES[jS])
                for w,pp in d.items():
                    q_=acc.get(w)
                    acc[w]=padd(q_,pp) if q_ is not None else pp
        Est[j]={w:pp for w,pp in acc.items() if pp}
    # Inn[W] = sum_c r_{c-1} prepend(c, Est[W-c])
    Inn=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for c in range(3, W+1):
            if (c-1) not in r: continue
            d=prepend(c, Est[W-c])
            for w,pp in d.items():
                q_=acc.get(w)
                pq=pmul(r[c-1], pp)
                acc[w]=padd(q_,pq) if q_ is not None else pq
        Inn[W]={w:pp for w,pp in acc.items() if pp}
    # Gg: log G~ coeffs: m L_m = (-1)^m[(a1+a3)^m+(a2+a4)^m-S4a^m]
    A13=[(0,1),(2,1)]
    Gg=[{} for _ in range(Wmax+1)]; Gg[0]={(): {ZERO:F(1)}}
    for j in range(1,Wmax+1):
        acc={}
        for m in range(2,j+1):
            lm=padd(linpow5(A13,m), linpow5(X24,m))
            lm=padd(lm, linpow5(S4,m), F(-1))
            lm=pscale(lm, F((-1)**m))
            d=st_letter(m, Gg[j-m])
            for w,pp in d.items():
                q_=acc.get(w)
                pq=pscale(pmul(pp, lm), F(1,j))
                acc[w]=padd(q_,pq) if q_ is not None else pq
        Gg[j]={w:pp for w,pp in acc.items() if pp}
    # A_reg[W] = sum st_dicts(Gg[i], Inn[j])
    Areg=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for i in range(0,W+1):
            j=W-i
            if not Gg[i] or not Inn[j]: continue
            d=st_dicts(Gg[i], Inn[j])
            for w,pp in d.items():
                q_=acc.get(w)
                acc[w]=padd(q_,pp) if q_ is not None else pp
        Areg[W]={w:pp for w,pp in acc.items() if pp}
    return {'E5':E5,'E24':E24,'ES':ES,'Est':Est,'Inn':Inn,'Gg':Gg,'Areg':Areg}

def rev_poly(p):
    return {(e[4],e[3],e[2],e[1],e[0]): c for e,c in p.items()}

def volreg_words(eng, Wmax=10):
    out=[{} for _ in range(Wmax+1)]
    for W in range(3,Wmax+1):
        acc={}
        for w,pp in eng['Areg'][W].items():
            acc[w]=padd(pp, rev_poly(pp))
        out[W]=acc
    return out

_TABLES={}
def table(W):
    if W not in _TABLES:
        import mzv_exact as M
        _TABLES[W]=M.DSTableFlint(W) if W>=11 else M.DSTable(W)
    return _TABLES[W]

def reduced_forms(Wmax=10):
    eng = build_engine(Wmax)
    vw = volreg_words(eng, Wmax)
    forms={}
    for W in range(3,Wmax+1):
        T=table(W); d=T.nullity
        acc=[{} for _ in range(d)]
        for w,pp in vw[W].items():
            assert sum(w)==W and w[0]>=2, (W,w)
            rw=T.reduce_word(w)
            for b in range(d):
                if rw[b]: acc[b]=padd(acc[b], pscale(pp, rw[b]))
        forms[W]=acc
    return forms, eng
'''
with open('engine5.py','w') as f:
    f.write(code)
print("engine5.py written", len(code))
