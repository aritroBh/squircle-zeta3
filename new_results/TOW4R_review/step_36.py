
import os, sys, pickle, time, math
import sympy as sp
from fractions import Fraction as Fr
base = "<project>"
work = base + "/research/08_theory/TOW4R_review"
sys.path.insert(0, base + "/research/08_theory/TOW2_tower_all_depths")
sys.path.insert(0, work + "/tow4_copy")
from fastred import QSI, IQ, AQ, BQ, ZERO
from genintegrand import S_poly
import towlib, tow4lib as T4

ONE = QSI(1)
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
    Binv0=B[0].inv(); out=[ZERO]*(N+1)
    for n in range(N+1):
        acc = A[n] if n<len(A) else ZERO
        for j in range(1,n+1):
            if j < len(B): acc = acc - B[j]*out[n-j]
        out[n]=acc*Binv0
    return out
def s_log(f, N):
    fp=[f[j+1].scale(Fr(j+1)) for j in range(len(f)-1)]
    q=s_div(fp,f,N-1)
    return [ZERO]+[q[r-1].scale(Fr(1,r)) for r in range(1,N+1)]
def s_pow_int(f,e):
    out=[ONE]
    for _ in range(e): out=s_mul(out,f)
    return out
N=8
tstar = QSI(Fr(1,2),0,Fr(1,2),0)
tser=[tstar, ZERO]+[ (IQ**(n+1)).scale(Fr(1,4)) for n in range(2,N+1)]
def logS_series(d, m, N):
    P = S_poly(d, m)
    S_ser=[ZERO]*(N+1)
    for s,c in enumerate(P):
        if c==0: continue
        tp=s_pow_int(tser,s)
        S_ser=s_add(S_ser, s_scale(tp+[ZERO]*(N+1-len(tp)), Fr(c)))
    return s_log(S_ser[:N+1], N)
def myL(d, m):
    L=s_log_pad = logS_series(d,m,N)
    L=[x for x in L]
    for r in range(1,N+1):
        L[r]=L[r]+(IQ**(-r)).scale(Fr((-1)**(r-1)*(2*d-1),r))
        L[r]=L[r]+((2*IQ)**(-r)).scale(Fr((-1)**(r-1)*(2-2*d),r))
        for pole in [AQ,-AQ,BQ,-BQ]:
            L[r]=L[r]+((IQ-pole)**(-r)).scale(Fr(-(-1)**(r-1),r))
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
            term=term.scale(Fr(1,math.factorial(mult)))
        tot=tot+term
    return tot

# exact Fraction evaluator for T4.arho
d_,m_=T4.d,T4.m
def make_eval(expr):
    num, den = sp.fraction(sp.together(expr))
    def prep(pe):
        P = sp.Poly(sp.expand(pe), d_, m_)
        return [(mon, Fr(sp.numer(sp.re(co))), Fr(sp.numer(sp.im(co)))) for mon, co in P.terms()]
    tn, td = prep(num), prep(den)
    def ev(terms, dv, mv):
        R=Fr(0); I=Fr(0)
        for mon, cr, ci in terms:
            f = Fr(dv)**mon[0]*Fr(mv)**mon[1]
            R += cr*f; I += ci*f
        return R,I
    def f(dv, mv):
        nr,ni = ev(tn,dv,mv); dr,di = ev(td,dv,mv)
        dd = dr*dr+di*di
        return ((nr*dr+ni*di)/dd, (ni*dr-nr*di)/dd)
    return f
ev7 = make_eval(T4.arho(7)); ev8 = make_eval(T4.arho(8))
t0=time.time()
badBell=0; badData=0; n=0
for dd in range(6,19):
    for p in range(dd,(dd-1)//2,-1):
        q=dd-p
        Lm=myL(dd,p)
        for delta, ev in [(7,ev7),(8,ev8)]:
            am = bell(Lm, delta)
            if am.p1!=0 or am.p3!=0: badBell+=1
            rr = ev(dd,p)
            if am.p0!=rr[0] or am.p2!=rr[1]:
                badBell+=1
                if badBell<4: print("Bell/T4 mismatch",(p,q),delta)
            # data ratio via towlib.rho (binteg)
            r0 = towlib.rho((p,q),0); rd = towlib.rho((p,q),delta)
            ratio = rd*r0.inv()
            if ratio.p1!=0 or ratio.p3!=0 or ratio.p0!=am.p0 or ratio.p2!=am.p2:
                badData+=1
                if badData<4: print("data mismatch",(p,q),delta, ratio, am)
        n+=1
print(f"a_7,a_8: myBell==T4.arho AND myBell==towlib.rho data ratio, {n} keys x2, badBell={badBell} badData={badData} ({round(time.time()-t0,1)}s)")
with open(work+"/progress.log","a") as fh:
    fh.write(f"{time.strftime('%H:%M:%S')} a_7/a_8 from-scratch Bell == T4.arho == towlib.rho(binteg) at {n} keys each; badBell={badBell} badData={badData}\n")
