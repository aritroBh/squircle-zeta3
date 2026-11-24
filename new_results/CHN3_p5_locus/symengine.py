
"""Symbolic engine: c_W(a) coefficient forms as dict-polys in a1..a4 (total degree W),
for each basis monomial. Polys: dict (e1,e2,e3,e4)->Fraction, homogeneous."""
from fractions import Fraction as F
import mzv_exact as M
import engine as E

def padd(p,q,s=F(1)):
    out=dict(p)
    for e,c in q.items():
        out[e]=out.get(e,F(0))+s*c
        if out[e]==0: del out[e]
    return out
def pscale(p,s):
    if s==0: return {}
    return {e:c*s for e,c in p.items()}
def pmul(p,q):
    out={}
    for e1,c1 in p.items():
        a1,a2,a3,a4=e1
        for e2,c2 in q.items():
            e=(a1+e2[0],a2+e2[1],a3+e2[2],a4+e2[3])
            out[e]=out.get(e,F(0))+c1*c2
    return out
def var(i,k):  # a_i^k
    e=[0,0,0,0]; e[i]=k; return {tuple(e):F(1)}
def linpow(terms, m):
    # (sum c_v a_v)^m, terms: list of (varindex, coeff)
    out={ (0,0,0,0):F(1) }
    # multinomial: iterate
    from math import comb, factorial
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
    rec(0,m,(0,0,0,0),F(1))
    return out

def symbolic_forms(Wmax=11):
    # base polys
    d=[{}]*(Wmax+2); g=[{}]*(Wmax+2); L=[{}]*(Wmax+2)
    for m in range(1,Wmax+2):
        sgn = 1 if m%2==1 else -1
        pA=linpow([(0,1),(2,1)],m)          # (a1+a3)^m
        pB=linpow([(0,1),(1,1),(2,1)],m)    # (a1+a2+a3)^m
        dm=padd(pA,pB,F(-1)); dm=pscale(dm,F(sgn,m))
        d[m]=dm
        g[m]=padd(dm, var(3,m), F(-1,m))
        sgn2 = 1 if m%2==0 else -1
        Lm=padd(var(1,m), pA); Lm=padd(Lm, pB, F(-1)); L[m]=pscale(Lm, F(sgn2,m))
    # delta
    delta=[{} for _ in range(Wmax+1)]; delta[0]={(0,0,0,0):F(1)}
    for l in range(1,Wmax+1):
        acc={}
        for i in range(1,l+1):
            acc=padd(acc, pscale(pmul(d[i],delta[l-i]), F(i,l)))
        delta[l]=acc
    # Fw[j]: word -> poly(deg j)
    Fw=[{} for _ in range(Wmax+1)]; Fw[0]={(): {(0,0,0,0):F(1)}}
    for j in range(1,Wmax+1):
        acc={}
        for i in range(1,j+1):
            for w,p in Fw[j-i].items():
                pg=pscale(pmul(p,g[i]), F(i,j))
                for w2,c2 in M.stuffle((i,),w).items():
                    q=acc.get(w2)
                    acc[w2] = padd(q,pg,c2) if q is not None else pscale(pg,c2)
        Fw[j]={w:p for w,p in acc.items() if p}
    # c0delta[s] = sum_{r+l=s-2} (-1)^{r+1} a3^{r+1} a4 * delta[l]
    c0d=[{} for _ in range(Wmax+1)]
    for s in range(2,Wmax+1):
        acc={}
        for r in range(0,s-1):
            l=s-2-r
            base=pscale(pmul(var(2,r+1),var(3,1)), F((-1)**(r+1)))
            acc=padd(acc, pmul(base,delta[l]))
        c0d[s]=acc
    # H[W]: word -> poly
    H=[{} for _ in range(Wmax+1)]; H[0]={(): {(0,0,0,0):F(1)}}
    for W in range(2,Wmax+1):
        acc={}
        for s in range(2,W+1):
            j=W-s
            for w,p in Fw[j].items():
                pq=pmul(c0d[s],p)
                acc[(s,)+w]=padd(acc.get((s,)+w,{}),pq)
        H[W]={w:p for w,p in acc.items() if p}
    # G[W]: word -> poly
    G=[{} for _ in range(Wmax+1)]; G[0]={(): {(0,0,0,0):F(1)}}
    for W in range(2,Wmax+1):
        acc={}
        for i in range(2,W+1):
            for w,p in G[W-i].items():
                pg=pscale(pmul(p,L[i]), F(i,W))
                for w2,c2 in M.stuffle((i,),w).items():
                    acc[w2]=padd(acc.get(w2,{}),pg,c2)
        G[W]={w:p for w,p in acc.items() if p}
    # final: forms[W][bindex] = poly, with reduction in the pair loop
    forms={}
    for W in range(2,Wmax+1):
        T=E.table(W); dnum=T.nullity
        Facc=[{} for _ in range(dnum)]
        for i in range(0,W+1):
            for w1,p1 in G[i].items():
                for w2,p2 in H[W-i].items():
                    st=M.stuffle(w1,w2)
                    # reduce the word-combination st
                    coords=[F(0)]*dnum
                    for w,c in st.items():
                        rw=T.reduce_word(w)
                        for b in range(dnum): coords[b]+=c*rw[b]
                    pq=pmul(p1,p2)
                    for b in range(dnum):
                        if coords[b]: Facc[b]=padd(Facc[b],pq,coords[b])
        forms[W]=Facc
    return forms

def peval(p, a):
    return sum(c*a[0]**e[0]*a[1]**e[1]*a[2]**e[2]*a[3]**e[3] for e,c in p.items())
