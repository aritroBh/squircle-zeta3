import sympy as sp, mpmath, pickle
from mpmath import mp, mpf, sqrt, zeta, log, nstr, factorial
mp.dps = 90
L2=log(2)
def eta(n): return (1-mpf(2)**(1-n))*zeta(n)
ck = pickle.load(open('ckpt_integrals90.pkl','rb'))
I  = {eval(k): mpf(v) for k,v in ck['I'].items()}
K2 = {eval(k): mpf(v) for k,v in ck['K2'].items()}
K3 = {eval(k): mpf(v) for k,v in ck['K3'].items()}
def K1(m): return 4*factorial(m)/mpf(2)**(m+1)*(eta(m+1)-eta(m))
def A5(a): return mpf(a+3)/9*factorial(a+5)*eta(a+6)/mpf(2)**(a+4)
def A7(a): return -factorial(a+7)*eta(a+8)/mpf(2)**(a+8)*(mpf(8*a)/81*(a+8)+mpf(488)/315)
R = sp.Rational
def Mval(a,b,c,d,e):
    if d==0 and e==0:
        if c==0: return I[(a,b)]
        if c==1: return -mpf(R(2*(a+3),3*(b+1)))*I[(a+2,b+1)]
        if c==2:
            return mpf(R(4,9))*(K1(a+6) if b==0 else (K2[a+6] if b==1 else K3[a+6]))
    if c==0 and e==0 and d==1:
        if b==0: return A5(a)
        m=a+6
        return (-mpf(R(2*(a+5),3*(b+1)))*I[(a+4,b+1)]
                + mpf(R(2,9))*(mpf(R(m*(m-1),b+1))*I[(m-2,b+1)] - b*(K1(m) if b==1 else K2[m])))
    if c==1 and d==1 and e==0 and b==0:
        return -mpf(R(2*(a+3),27))*K1(a+8)
    if e==1 and b==0 and c==0 and d==0:
        return A7(a)
    if e==1 and b==1 and c==0 and d==0:
        m1,m2,m3 = a+9,a+8,a+7
        return (mpf(R(4,81))*(-mpf(R(m1*(m1-1)*(m1-2),2))*I[(m1-3,2)] + mpf(R(3*m1,2))*K1(m1-1))
                + mpf(R(4,9))*(mpf(R(m2*(m2-1),2))*I[(m2-2,2)] - K1(m2))
                + mpf(R(244,315))*(-mpf(R(m3,2))*I[(m3-1,2)]))
    raise ValueError((a,b,c,d,e))
mons = pickle.load(open('ckpt_monomials.pkl','rb'))
F = {}
for k in range(1,9):
    tot = mpf(0)
    for (cs,a,b,c,d,e) in mons[k]:
        tot += mpf(cs)*Mval(a,b,c,d,e)
    F[k] = tot
    print(f"F_{k} = {nstr(F[k],45)}   pred c{k+1} = {nstr(4*sqrt(2)*F[k],45)}")
pickle.dump({k: nstr(v,80) for k,v in F.items()}, open('ckpt_Fk_assembled.pkl','wb'))
