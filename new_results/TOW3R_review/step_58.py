work = '<project>/research/08_theory/TOW3R_review/'
src = '''import sys, os, time, pickle
from fractions import Fraction as Fr
from math import comb, factorial
work = os.path.dirname(os.path.abspath(__file__))

# ---------- exact Gaussian-rational + QSI arithmetic ----------
# Gaussian rational: (re, im) of Fractions. QSI: (A, B) with value A + B*sqrt2, A,B Gaussian.
I2 = (Fr(0), Fr(1))
def gadd(x,y): return (x[0]+y[0], x[1]+y[1])
def gsub(x,y): return (x[0]-y[0], x[1]-y[1])
def gmul(x,y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def gdiv(x,y):
    den = y[0]**2+y[1]**2
    return ((x[0]*y[0]+x[1]*y[1])/den, (x[1]*y[0]-x[0]*y[1])/den)
def gpow_i(n):  # i^n
    return [(Fr(1),Fr(0)),(Fr(0),Fr(1)),(Fr(-1),Fr(0)),(Fr(0),Fr(-1))][n%4]
def qadd(x,y): return (gadd(x[0],y[0]), gadd(x[1],y[1]))
def qsub(x,y): return (gsub(x[0],y[0]), gsub(x[1],y[1]))
def qmul(x,y):  # (A+Bs2)(C+Ds2) = AC+2BD + (AD+BC)s2
    return (gadd(gmul(x[0],y[0]), gmul(gmul(x[1],y[1]),(Fr(2),Fr(0)))), gadd(gmul(x[0],y[1]), gmul(x[1],y[0])))
def qscale(x, c): return (gmul(x[0],c), gmul(x[1],c))  # c Gaussian
def qpow(x, n):
    r = ((Fr(1),Fr(0)),(Fr(0),Fr(0)))
    for _ in range(n): r = qmul(r, x)
    return r
def qconj(x): return ((x[0][0],-x[0][1]),(x[1][0],-x[1][1]))

S2 = (Fr(2),Fr(0))
aq = ((Fr(1),Fr(0)),(Fr(1),Fr(0)))   # a = 1 + sqrt2
bq = ((Fr(-1),Fr(0)),(Fr(1),Fr(0)))  # b = -1 + sqrt2
def tau(n, which):
    # s = n+1; (-1)^s [ (i - o1)^{-(s+1)} + (i + o2)^{-(s+1)} - i^{-(s+1)} ]
    s = n+1
    iQ = ((Fr(0),Fr(1)),(Fr(0),Fr(0)))
    if which==1: o1, o2 = bq, aq   # (i-b), (i+a)
    else: o1, o2 = aq, bq          # (i-a), (i+b)
    t1 = qpow(( (Fr(0),Fr(1)), (Fr(-1)*o1[1][0], Fr(0)) ), 0) # placeholder
    # compute (i - o1)^{-(s+1)}: i-o1 as QSI
    def qinv(x):
        # 1/(A+Bs2) = (A - B s2)/(A^2 - 2 B^2)
        A,B = x
        den = gsub(gmul(A,A), gscale(gmul(B,B),(Fr(2),Fr(0))))
        return (gdiv(A,den), gdiv((-B[0],-B[1]),den))
    imo1 = ((Fr(0)-o1[0][0], Fr(1)-o1[0][1]), (Fr(-1)*o1[1][0], Fr(-1)*o1[1][1]))
    ipo2 = ((Fr(0)+o2[0][0], Fr(1)+o2[0][1]), o2[1])
    r = qadd(qpow(qinv(imo1), s+1), qpow(qinv(ipo2), s+1))
    # subtract i^{-(s+1)} as QSI
    iv = gpow_i(-(s+1))
    r = qsub(r, (iv,(Fr(0),Fr(0))))
    if s%2==1: r = qscale(r, (Fr(-1),Fr(0)))
    return r
_tau = {}
def tau1(n):
    if ('1',n) not in _tau: _tau[('1',n)] = tau(n,1)
    return _tau[('1',n)]
def tau2(n):
    if ('2',n) not in _tau: _tau[('2',n)] = tau(n,2)
    return _tau[('2',n)]

# sanity vs tow2lib table
sys.path.insert(0, work+'/tow3_copy')
import sympy as sp, tow2lib as T
for n in range(8):
    for w, f in [(1,tau1),(2,tau2)]:
        v = f(n)
        assert v[1]==(Fr(0),Fr(0)), (n,w,v)  # tau in Q(i)
        tv = T.tau1f(n) if w==1 else T.tau2f(n)
        assert (v[0][0], v[0][1]) == (sp.Rational(sp.re(tv)), sp.Rational(sp.im(tv))), (n,w,v,tv)
print('tau1/tau2 from definition == tow2lib table (n=0..7), Q(i)-valued: OK', flush=True)

# ---------- tower data ----------
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
def u(p,q,k,dl):
    x = H[(p,q)][dl][k]
    return ((x[0],x[2]),(x[1],x[3]))  # ((Ar,Ai),(Br,Bi))

def fold(p, delta, k):
    # Fold_p(delta; u_{pp}^{(k)}) : -p sum_{mu=m}^{4p-4} (-i)^mu binom(mu-2, mu-m) i^{m-2} conj(g_mu), m=4p-4-delta
    m = 4*p-4-delta
    tot = ((Fr(0),Fr(0)),(Fr(0),Fr(0)))
    for mu in range(m, 4*p-3):
        t = 4*p-4-mu
        g = ((Fr(0),Fr(0)),(Fr(0),Fr(0)))
        for dp in range(0, t+1):
            uu = u(p,p,k,dp)
            c = qscale(tau1(t-dp), (Fr(1)/ (4*p-3-dp), Fr(0)))
            g = qsub(g, qmul(uu, c))
        gc = qconj(g)
        coef1 = gmul(gpow_i(mu), (Fr(-1)**mu if mu%2 else Fr(1), Fr(0))) # (-i)^mu
        coef1 = gmul(coef1, (Fr(comb(mu-2, mu-m)), Fr(0)))
        coef1 = gmul(coef1, gpow_i(m-2))
        tot = qadd(tot, qscale(gc, coef1))
    return qscale(tot, (Fr(-p),Fr(0)))

def T1_rhs(p, q, k, delta):
    d = p+q
    tot = ((Fr(0),Fr(0)),(Fr(0),Fr(0)))
    for dp in range(0, delta+1):
        c1 = qscale(tau1(delta-dp), (Fr(p+1, 4*0+1)*0+Fr(p+1,1)/(2*d-dp-1), Fr(0)))
        tot = qadd(tot, qmul(u(p+1,q,k,dp), c1))
        if p > q:
            c2 = qscale(tau2(delta-dp), (Fr(q+1,1)/(2*d-dp-1), Fr(0)))
            tot = qadd(tot, qmul(u(p,q+1,k,dp), c2))
    if q == p-1:
        tot = qadd(tot, fold(p, delta, k))
    return tot

t0=time.time()
bad=[]; tested=0
for p in range(1, 40):
    for q in range(0, p+1):
        d = p+q
        if d < 12 or d > 30: continue
        if (p,q) not in H: continue
        dmax = min(6, 2*d-4)
        for k in range(d+1, 53):
            for dl in range(0, dmax+1):
                # need parent data
                try:
                    rhs = T1_rhs(p,q,k,dl)
                except KeyError:
                    continue
                lhs = u(p,q,k,dl)
                tested += 1
                if lhs != rhs:
                    bad.append(((p,q),k,dl,lhs,rhs))
print(f'T1 recursion at fresh cells (degrees 12..30, k=d+1..52, depths<=6): tested {tested}, bad {len(bad)}  ({time.time()-t0:.0f}s)', flush=True)
for b in bad[:5]: print(b)
'''
open(work+'step_14.py','w').write(src)
import subprocess, datetime
PY='<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,work+'step_14.py'],capture_output=True,text=True,timeout=1700)
print(r.stdout, r.stderr[-2500:])
open(work+'step_14.out','w').write(r.stdout+r.stderr)
with open(work+'progress.log','a') as f:
    f.write(f"{datetime.datetime.now().isoformat()} step_14 own exact T1 recursion vs fresh data degrees 12..30\n")
