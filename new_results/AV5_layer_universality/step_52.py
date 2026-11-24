import sympy as sp, mpmath, pickle
from mpmath import mp, mpf, sqrt, zeta, log, polylog, nstr, factorial
mp.dps = 200
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
ckA = pickle.load(open('ckpt_integrals135.pkl','rb'))
ckB = pickle.load(open('ckpt_integrals210.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ckA['I'].items()}
for k,v in ckB['I'].items(): I[eval(k)] = mpf(v)
K2 = {eval(k): mpf(v) for k,v in ckA['K2'].items()}
for k,v in ckB['K2'].items(): K2[eval(k)] = mpf(v)
K3 = {eval(k): mpf(v) for k,v in ckA['K3'].items()}
def K1(m): return 4*factorial(m)/mpf(2)**(m+1)*(eta(m+1)-eta(m))
def A5(a): return mpf(a+3)/9*factorial(a+5)*eta(a+6)/mpf(2)**(a+4)
def A7(a): return -factorial(a+7)*eta(a+8)/mpf(2)**(a+8)*(mpf(8*a)/81*(a+8)+mpf(488)/315)
R = sp.Rational
def Iget(a,b):
    if b==1: return factorial(a)*eta(a+2)/mpf(2)**(a+1)
    return I[(a,b)]
def Mval(a,b,c,d,e):
    if d==0 and e==0:
        if c==0: return Iget(a,b)
        if c==1: return -mpf(R(2*(a+3),3*(b+1)))*Iget(a+2,b+1)
        if c==2: return mpf(R(4,9))*(K1(a+6) if b==0 else (K2[a+6] if b==1 else K3[a+6]))
    if c==0 and e==0 and d==1:
        if b==0: return A5(a)
        m=a+6
        return (-mpf(R(2*(a+5),3*(b+1)))*Iget(a+4,b+1)
                + mpf(R(2,9))*(mpf(R(m*(m-1),b+1))*Iget(m-2,b+1) - b*(K1(m) if b==1 else K2[m])))
    if c==1 and d==1 and e==0 and b==0: return -mpf(R(2*(a+3),27))*K1(a+8)
    if e==1 and b==0 and c==0 and d==0: return A7(a)
    if e==1 and b==1 and c==0 and d==0:
        m1,m2,m3 = a+9,a+8,a+7
        return (mpf(R(4,81))*(-mpf(R(m1*(m1-1)*(m1-2),2))*Iget(m1-3,2) + mpf(R(3*m1,2))*K1(m1-1))
                + mpf(R(4,9))*(mpf(R(m2*(m2-1),2))*Iget(m2-2,2) - K1(m2))
                + mpf(R(244,315))*(-mpf(R(m3,2))*Iget(m3-1,2)))
    raise ValueError((a,b,c,d,e))
mons = pickle.load(open('ckpt_monomials.pkl','rb'))
F = {}
for k in [7,8]:
    tot = mpf(0)
    for (cs,a,b,c,d,e) in mons[k]:
        tot += mpf(cs)*Mval(a,b,c,d,e)
    F[k] = tot
    print(f"F_{k} = {nstr(F[k],60)}")
pickle.dump({k: nstr(v,190) for k,v in F.items()}, open('ckpt_Fk78_200.pkl','wb'))
# zeta-only PSLQ at dps 200
Z = lambda n: zeta(n)
B8 = [('Z2^4',Z(2)**4),('Z2Z3^2',Z(2)*Z(3)**2),('Z3Z5',Z(3)*Z(5)),('Z8',Z(8))]
B9 = [('Z2^3Z3',Z(2)**3*Z(3)),('Z2^2Z5',Z(2)**2*Z(5)),('Z2Z7',Z(2)*Z(7)),('Z3^3',Z(3)**3),('Z9',Z(9))]
for k,B in [(7,B8),(8,B9)]:
    rel = mpmath.pslq([v for _,v in B]+[F[k]], maxcoeff=10**12, maxsteps=50000)
    print(f"F_{k} zeta-only PSLQ (dps200, maxcoeff 1e12):", rel)
