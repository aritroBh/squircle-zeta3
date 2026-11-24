import pathlib
code = r'''
import sys, time, math, pickle
sys.dont_write_bytecode = True
import mpmath as mp
from fractions import Fraction as Fr
import sympy as sp

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t0 = time.time()
log("rul1_step04b: N4, N5 from D12's exact formula (PROVED) + proved J-values")

# ---- rebuild PSI[0..4] exactly (D12 step_04 construction, rational arithmetic) ----
KMAX = 5
def addC(A,B):
    C=dict(A)
    for k,v in B.items():
        C[k]=C.get(k,Fr(0))+v
        if C[k]==0: del C[k]
    return C
def scaleC(A,c): return {k:v*c for k,v in A.items()}
def mulC(A,B):
    C={}
    for (i1,j1,k1),v1 in A.items():
        for (i2,j2,k2),v2 in B.items():
            key=(i1+i2,j1+j2,k1+k2); C[key]=C.get(key,Fr(0))+v1*v2
    return {k:v for k,v in C.items() if v}
def serMul(A,B,K=KMAX):
    C=[{} for _ in range(K+1)]
    for i in range(len(A)):
        if not A[i]: continue
        for j in range(len(B)):
            if not B[j]: continue
            if i+j>K: break
            C[i+j]=addC(C[i+j],mulC(A[i],B[j]))
    return C
def serExp(A,K=KMAX):
    R=[{} for _ in range(K+1)]; R[0]={(0,0,0):Fr(1)}
    for n in range(1,K+1):
        acc={}
        for k in range(1,n+1):
            if A[k]: acc=addC(acc,mulC(scaleC(A[k],Fr(k)),R[n-k]))
        R[n]=scaleC(acc,Fr(1,n))
    return R
y=sp.symbols('y')
lc=sp.log(sp.cosh(y)).series(y,0,KMAX+2).removeO()
c=[None]*(KMAX+1); c[1]=Fr(-1,2)
for r in range(2,KMAX+1): c[r]=Fr(sp.expand(lc).coeff(y,r))/Fr(2)**r
T=[{} for _ in range(KMAX+1)]; T[0]={(0,0,1):Fr(1)}
for r in range(1,KMAX+1):
    if c[r]==0: continue
    for t in range(0,KMAX+1-r):
        T[r+t]=addC(T[r+t],{(r,0,0):Fr(math.comb(r+t-1,t))*c[r]})
Q=[{} for _ in range(KMAX+1)]
for n in range(KMAX+1):
    Q[n]=addC(Q[n],scaleC(T[n],Fr(-1)))
    if n+1<=KMAX: Q[n+1]=addC(Q[n+1],scaleC(T[n],Fr(-1)))
psi2=[scaleC(co,Fr(1,2)) for co in serExp([{}]+Q[1:])]
B=[{} for _ in range(KMAX+1)]
for m in range(2,KMAX+1): B[m]={(1,0,0):Fr(-1)}
psi1=serExp(B)
psi3=[{} for _ in range(KMAX+1)]
for m in range(1,KMAX+1): psi3[m]={(0,m,0):Fr(1,math.factorial(m))}
PSI=serMul(serMul(psi1,psi2),psi3)
for m in range(KMAX+1):
    print("PSI[%d] = %s" % (m, PSI[m]))

# ---- J-values as exact named-basis data ----
import rul1, rul1_reduce, rul1num
certs = pickle.load(open('rul1_certs_J.pkl','rb'))
l2, pi2, z3, z5, z2 = sp.symbols('l2 pi2 z3 z5 z2')
# zeta values in terms of our symbols: zeta(2)=pi2/6, zeta(4)=pi2^2/90
def zeta(n):
    if n == 2: return pi2/6
    if n == 3: return z3
    if n == 4: return pi2**2/90
    if n == 5: return z5
    raise ValueError(n)
def eta(n):
    return (1-Fr(2)**(1-n))*zeta(n)
J = {}
for k in range(0, 10):
    J[(1,k)] = math.factorial(k)*eta(k+2)
J[(2,0)] = z3/4
for (a,k),(cc,h) in certs.items():
    expr = 0
    for nm, x in cc.items():
        expr += x*sp.symbols(nm)
    J[(a,k)] = expr
# symbols for named constants
namedsym = {nm: sp.symbols(nm) for nm in ['log2^4','pi2*log2^2','pi4','z3*log2','Li4h','log2^5','pi2*log2^3','z3*log2^2','z3*pi2','Li5h','Li4h*log2']}
sub_base = {'Li4h': sp.symbols('Li4h'), 'Li5h': sp.symbols('Li5h')}

def Jexpr(a,k):
    return J[(a,k)]

# ---- assemble INT[m], N4, N5 symbolically ----
INT = []
for m in range(KMAX+1):
    s = 0
    for (i,j,k), v in PSI[m].items():
        s += v*Jexpr(j,i)*l2**k
    INT.append(sp.expand(s))
N4 = sp.expand(Fr(-1,1)**5*l2**5/120 + sum(INT[0:4]))
N5 = sp.expand(Fr(-1,1)**6*l2**6/720 + sum(INT[0:5]))
print("N4 (raw) =", N4)
print("N5 (raw) =", N5)

# substitute pi4 = pi2^2/90 (our named pi4 symbol) -- the J-closed forms used pi4 symbol;
# express everything in the small symbol set
subs = {'pi4': pi2**2/90, 'z3*pi2': z3*pi2, 'log2^4': l2**4, 'log2^5': l2**5,
        'pi2*log2^2': pi2*l2**2, 'pi2*log2^3': pi2*l2**3, 'z3*log2': z3*l2,
        'z3*log2^2': z3*l2**2, 'Li4h': sp.symbols('Li4h'), 'Li5h': sp.symbols('Li5h')}
N4s = sp.expand(N4.subs(subs, simultaneous=True))
N5s = sp.expand(N5.subs(subs, simultaneous=True))
print("N4 =", N4s)
print("N5 =", N5s)

# expected closed forms (D12 result.md sec 1.3)
N4_exp = -l2**5/120 + pi2/24*(1-l2+l2**2/2) + z3/4*(1-l2) + pi2**2/(90*2880)*90
N4_exp = -l2**5/120 + pi2/24*(1-l2+l2**2/2) + z3/4*(1-l2) + pi2**2/(90*2880/90)
# careful: pi^4/2880 = (pi2^2)/2880
N4_exp = -l2**5/120 + pi2/24*(1-l2+l2**2/2) + z3/4*(1-l2) + pi2**2/2880
N5_exp = (l2**6/720 + pi2/24*(1-l2+l2**2/2-l2**3/6) + z3/4*(1-l2+l2**2/2)
          + pi2**2/2880*(1-l2) - 7*pi2**2/(1920*90)*90 - z3*pi2/48 + z5/8)
# careful: 7 pi^4/1920 = 7 pi2^2/1920
N5_exp = (l2**6/720 + pi2/24*(1-l2+l2**2/2-l2**3/6) + z3/4*(1-l2+l2**2/2)
          + pi2**2/2880*(1-l2) - 7*pi2**2/1920 - z3*pi2/48 + z5/8)
d4 = sp.expand(N4s - N4_exp)
d5 = sp.expand(N5s - N5_exp)
print("N4 - expected =", d4)
print("N5 - expected =", d5)
assert d4 == 0 and d5 == 0
log("TARGET 1c PASS: N4, N5 closed forms follow EXACTLY from proved formula + proved J-values")

# ---- numeric verification at 130 digits: direct fit of pi_p ----
mp.mp.dps = 150
def pi_h(h):
    f = lambda u: mp.cosh(u/2)**(-1-h)*mp.cosh(u*(1-h)/(2*h))**h
    return 2*h*mp.quad(f,[0,mp.inf])/4
ps=[70+4*i for i in range(12)]
fs=[pi_h(1/mp.mpf(p)) for p in ps]
A=mp.matrix(12,7); b=mp.matrix(12,1)
for a,p in enumerate(ps):
    h=1/mp.mpf(p)
    for k in range(7): A[a,k]=h**(k+1)
    b[a]=fs[a]-1
coef=mp.lu_solve(A,b)
NV = rul1num.named_values(150)
def evalN(Nexpr):
    return mp.mpf(str(Nexpr.subs({l2: NV['log2'], pi2: NV['pi2'], z3: NV['z3'], z5: NV['z5'],
                                  sp.symbols('Li4h'): NV['Li4h'], sp.symbols('Li5h'): NV['Li5h']})))
for k,(Nexpr) in [(4,N4s),(5,N5s)]:
    fit = coef[k]; exact = evalN(Nexpr)
    print("N_%d: direct-fit = %s" % (k, mp.nstr(fit, 25)))
    print("     closed form = %s" % mp.nstr(exact, 25))
    print("     |diff| = %s" % mp.nstr(abs(fit-exact),3))
log("TARGET 1d: closed forms vs direct pi_p fit printed (expect ~1e-100 agreement)")
log("rul1_step04b done %.1fs" % (time.time()-t0))
'''
pathlib.Path('rul1_step04b.py').write_text(code)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice',PY,'rul1_step04b.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout[-4500:]); print(r.stderr[-1200:])
