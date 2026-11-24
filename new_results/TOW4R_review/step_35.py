
import os, sys, pickle, time
import sympy as sp
from fractions import Fraction as Fr
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, base + "/research/08_theory/TOW2_tower_all_depths")
sys.path.insert(0, work + "/tow4_copy")
from fastred import QSI, IQ, AQ, BQ, ZERO
from genintegrand import S_poly
import towlib, tow2lib as T2, tow4lib as T4

ONE = QSI(1)
def ser_const(c): return [c]
def s_add(A,B):
    n = max(len(A),len(B)); out=[ZERO]*n
    for j in range(n):
        out[j] = (A[j] if j<len(A) else ZERO)+(B[j] if j<len(B) else ZERO)
    return out
def s_scale(A,q): return [a*q for a in A]
def s_mul(A,B):
    out=[ZERO]*(len(A)+len(B)-1)
    for a,ca in enumerate(A):
        if ca.is_zero(): continue
        for b,cb in enumerate(B):
            if not cb.is_zero(): out[a+b]=out[a+b]+ca*cb
    return out
def s_div(A,B,N):
    # A/B series to order N (B[0] != 0)
    Binv=[ZERO]*(N+1); Binv[0]=B[0].inv()
    out=[ZERO]*(N+1)
    for n in range(N+1):
        acc = A[n] if n<len(A) else ZERO
        for j in range(1,n+1):
            if j < len(B): acc = acc - B[j]*out[n-j]
        out[n]=acc*Binv[0]
    return out
def s_log(f, N):
    # [w^r] log f, r=1..N, given f series to order N, f[0]!=0 : via f'/f
    fp=[f[j+1].scale(Fr(j+1)) for j in range(len(f)-1)]
    q=s_div(fp,f,N-1)
    return [ZERO]+[q[r-1].scale(Fr(1,r)) for r in range(1,N+1)]
def s_pow_int(f,e):
    out=[ONE]
    for _ in range(e): out=s_mul(out,f)
    return out

N=8
# t(i+w) series: t* + sum_{n>=2} i^{n+1}/4 w^n   (t'(i)=0)
tstar = QSI(Fr(1,2),0,Fr(1,2),0)
tser=[tstar, ZERO]+[ (IQ**(n+1)).scale(Fr(1,4)) for n in range(2,N+1)]
w_ser=[ZERO,ONE]  # v = i + w
v_ser=[QSI(0,0,1,0),ONE]  # i + w

def logS_series(d, m, N):
    P = S_poly(d, m)   # rational coeffs low-to-high
    # S(t(i+w)): compose
    S_ser=[ZERO]*(N+1)
    for s,c in enumerate(P):
        if c==0: continue
        tp=s_pow_int(tser,s)
        S_ser=s_add(S_ser, s_scale(tp+[ZERO]*(N+1-len(tp)), Fr(c)))
    S_ser=S_ser[:N+1]
    return s_log(S_ser, N), S_ser

def myL(d, m):
    # Phi(v) = v^{2d-1} S(t(v)) (v+i)^{2-2d} / ((v^2-a^2)(v^2-b^2))
    L=[ZERO]*(N+1)
    lg,_ = logS_series(d,m,N)
    L=s_add(L,lg)
    # (2d-1) log v :  f = v,  f'/f = 1/v
    logv=[ZERO]*(N+1)
    for r in range(1,N+1):
        # d^r/dv^r log v /r! at i = (-1)^{r-1}/(r i^r)
        logv[r]=(IQ**(-r)).scale(Fr((-1)**(r-1),r))
    L=s_add(L,s_scale(logv,Fr(2*d-1)))
    # (2-2d) log(v+i): (-1)^{r-1}/(r (2i)^r)
    logvi=[ZERO]*(N+1)
    for r in range(1,N+1):
        logvi[r]=((2*IQ)**(-r)).scale(Fr((-1)**(r-1),r))
    L=s_add(L,s_scale(logvi,Fr(2-2*d)))
    # -log(v^2-a^2)-log(v^2-b^2): -(-1)^{r-1}/r * sum (v0)^{-r}
    for pole in [AQ,-AQ,BQ,-BQ]:
        lg=[ZERO]*(N+1)
        for r in range(1,N+1):
            lg[r]=((IQ-pole)**(-r)).scale(Fr((-1)**(r-1),r))
        L=s_add(L,s_scale(lg,Fr(-1)))
    return L

def bell(L, delta):
    from collections import Counter
    def parts(rem, mn, cur):
        if rem==0:
            yield tuple(sorted(cur)); return
        for p in range(mn, rem+1):
            yield from parts(rem-p,p,cur+[p])
    tot=ZERO
    for part in parts(delta,1,[]):
        C=Counter(part); term=ONE
        for sz,mult in C.items():
            for _ in range(mult): term=term*L[sz]
            term=term.scale(Fr(1,__import__('math').factorial(mult)))
        tot=tot+term
    return tot

d_,m_=T2.d,T2.m
t0=time.time()
keys=[]
for dd in range(6,19):
    for p in range(dd,(dd-1)//2,-1):
        keys.append((p,dd-p))
print("n keys:", len(keys))
badL=0; badA=0; badD=0; n=0
for (p,q) in keys:
    dd=p+q
    Lm=myL(dd,p)
    for r in range(1,9):
        ref=T2.L(r).subs({d_:dd,m_:p})
        # ref is sympy expr in Q(i,sqrt2)? tow2lib L uses aq,bq sqrt2 -> sympy simplify to QSI compare
        ref=sp.nsimplify(ref) if False else sp.simplify(ref)
        # convert sympy to (re,im) rationals if in Q(i); L may contain sqrt2? p5 sums cancel sqrt2.
        if ref.has(sp.sqrt(2)):
            print("sqrt2 remains in tow2lib.L!", r); break
        rr=(Fr(sp.numer(sp.re(ref)))/Fr(sp.denom(sp.re(ref))), Fr(sp.numer(sp.im(ref)))/Fr(sp.denom(sp.im(ref))))
        got=Lm[r]
        if got.p1!=0 or got.p3!=0 or got.p0!=rr[0] or got.p2!=rr[1]:
            badL+=1
            if badL<4: print("L mismatch", (p,q), r, got, rr)
    n+=1
print(f"L_r from-scratch vs tow2lib: {n} keys x 8, bad = {badL}  ({round(time.time()-t0,1)}s)")
import pickle
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} from-scratch L_r (genintegrand S_poly + QSI series) vs tow2lib.L: {n} keys x r=1..8, bad={badL}\n")
