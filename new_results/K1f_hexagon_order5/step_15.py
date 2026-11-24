
import sympy as sp, pickle, time
WD = '<project>/research/08_theory/K1f_hexagon_order5'

# ---- General single-ratio edge TPS engine (E1/E5 style: Fpow * sqrt(M)) ----
A0,A1,B,kap = sp.symbols('A0 A1 B kap')
L2,q0,V02 = sp.symbols('L2 q0 V02')
lA0,lA1,lB = sp.symbols('lA0 lA1 lB')
t,X,Y,W = sp.symbols('t X Y W')
N = 6

def mul(a,b):
    return [sp.expand(sum(a[i]*b[k-i] for i in range(k+1))) for k in range(N)]
def inv(a):
    r=[sp.Integer(1)]+[0]*(N-1)
    for k in range(1,N):
        r[k]=sp.expand(-sum(a[i]*r[k-i] for i in range(1,k+1)))
    return r
def log1(a):
    x=[sp.Integer(0)]+a[1:]
    res=[sp.Integer(0)]*N; pw=[sp.Integer(1)]+[0]*(N-1)
    for j in range(1,N):
        pw=mul(pw,x); res=[sp.expand(res[k]+sp.Rational((-1)**(j+1),j)*pw[k]) for k in range(N)]
    return res
def exp0(a):
    res=[sp.Integer(1)]+[0]*(N-1); pw=[sp.Integer(1)]+[0]*(N-1)
    for j in range(1,N):
        pw=mul(pw,a); res=[sp.expand(res[k]+pw[k]/sp.factorial(j)) for k in range(N)]
    return res

a0=A0*t; a1=A1*(1-t); a2=B*(1+kap*t)
LA0=X+lA0; LA1=Y+lA1; LA2=W+lB
t0=time.time()
T=[sp.Integer(1)]+[sp.expand(a0*LA0**k+a1*LA1**k+a2*LA2**k) for k in range(1,N)]
M=[sp.Integer(0)]+[sp.expand(A0*LA0**k-A1*LA1**k+B*kap*LA2**k) for k in range(1,N)]
F=[sp.expand(T[k]/sp.factorial(k)) for k in range(N)]
Num=[sp.expand(M[k]/sp.factorial(k)) for k in range(N)]
Phi=mul(Num,inv(F))
oneoverp=[sp.Integer((-1)**k) for k in range(N)]
Fpow=exp0([sp.expand(-c) for c in mul(oneoverp,log1(F))])
# metric
yyp=q0+L2*t; yy=V02+2*q0*t+L2*t**2
Phi2=mul(Phi,Phi)
Mmet=[sp.Integer(0)]*N
# M(d) = L2 - 2 yyp Phi + yy Phi^2 ; Phi[0]=0
for k in range(N):
    val = -2*yyp*Phi[k] + yy*Phi2[k]
    if k==0: val += L2
    Mmet[k]=sp.expand(val)
# sqrt(M) = sqrt(L2)*sqrt(1+Xmetric), Xmetric=(M-L2)/L2
Xmet=[sp.Integer(0)]+[sp.expand(Mmet[k]/L2) for k in range(1,N)]
sq=[sp.Integer(1)]+[0]*(N-1); pw=[sp.Integer(1)]+[0]*(N-1)
for j in range(1,N):
    pw=mul(pw,Xmet); sq=[sp.expand(sq[k]+sp.binomial(sp.Rational(1,2),j)*pw[k]) for k in range(N)]
Smet=[sp.expand(sp.sqrt(L2)*c) for c in sq]
E_a=mul(Fpow,Smet); E_a=[sp.expand(c) for c in E_a]
print("Approach (a) built, time %.1f s"%(time.time()-t0))
print("E5 (a) terms:", len(sp.Add.make_args(E_a[5])))
pickle.dump({'E':[str(c) for c in E_a]}, open(f'{WD}/E5_general_approachA.pkl','wb'))
print("saved E5_general_approachA.pkl")
