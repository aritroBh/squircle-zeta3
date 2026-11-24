
import mpmath as mp
mp.mp.dps = 80
def A_val(u1,u2,u3,u4,u5,u,M=4000):
    # A(u1,u2;u4,u5) at direction*u, with tail Richardson
    S4 = (u1+u2+u3+u4)*u
    def alpha(m):
        return mp.rf(-u5*u,m)/mp.factorial(m)*(u4*u)/(u4*u+m)*u2*u*mp.beta(u2*u+u4*u+m, 1+u1*u+u3*u)
    Ms=[M, M*3//4, M//2, M*5//12, M//3]
    vals=[mp.fsum(alpha(m) for m in range(Mi)) for Mi in Ms]
    s_exp = (u1+u3+u5)*u+3
    n=len(Ms)
    Am=mp.matrix(n,n)
    for i,Mi in enumerate(Ms):
        Am[i,0]=1
        for j in range(1,n): Am[i,j]=Mi**(1-s_exp)*Mi**(-(j-1))
    return mp.lu_solve(Am, mp.matrix(vals))[0]

# A at (0,1,1,1,1)*u, extract c2..c5
hs=[mp.mpf(1)/200, mp.mpf(1)/400, mp.mpf(1)/800, mp.mpf(1)/1600, mp.mpf(1)/3200, mp.mpf(1)/6400]
Vs=[A_val(0,1,1,1,1,h) for h in hs]
n=len(hs)
Amat=mp.matrix(n,n)
for i,h in enumerate(hs):
    for j in range(n): Amat[i,j]=h**(j+1)   # A = a2/(a2+a4) + c2 u^2... wait m0 ~ 1/2; include u^1? A has no u^1
b=mp.matrix([V-mp.mpf('0.5') for V in Vs])
# fit c2..c6 (5 unknowns) using 6 points least squares? use 5 of them
Amat=mp.matrix(5,5)
for i in range(5):
    for j in range(5): Amat[i,j]=hs[i]**(j+2)
sol=mp.lu_solve(Amat, mp.matrix([Vs[i]-mp.mpf('0.5') for i in range(5)]))
print("A c2..c6:", [mp.nstr(s,14) for s in sol])
print("A c5/z5:", mp.nstr(sol[3]/mp.zeta(5),14))
# m0A c5 = -21 zeta(5) -> regular c5 = A_c5 +21 z5
print("regular c5 (G~Sa) = A_c5 + 21*z5 =", mp.nstr(sol[3]+21*mp.zeta(5),16))
print("true regular c5 = -15*z5 + 2*z2*z3 =", mp.nstr(-15*mp.zeta(5)+2*mp.zeta(2)*mp.zeta(3),16))
