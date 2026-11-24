code = r'''
from squircle import *
from mpmath import mp, zeta, mpf, nstr, fabs, fsum, power, dirichlet
setdps(40)
p = mpf(2); s = mpf(4); a = s/p
def h(m,n): return (m**p + n**p)**(-a)
N, M = 30, 300
def T1brute(m, K2=4000):
    return fsum([h(m,n) for n in range(1,K2+1)]) + power(m,1-s)*Jint(mpf(K2)/m,p,a)
E1_brute = fsum([T1brute(m) for m in range(N+1, M+1)])
tailM = Aint(p,a)*harmtail(M, s-1) - harmtail(M, s)/2
E1_brute += tailM
exact = 4*zeta(2)*dirichlet(2,[0,1,0,-1])
S_exact = (exact - 4*zeta(s))/4
Sc = fsum([h(m,n) for m in range(1,N+1) for n in range(1,N+1)])
print("S exact :", nstr(S_exact, 25))
print("Sc+2E1  :", nstr(Sc + 2*E1_brute, 25))
print("diff    :", nstr(Sc + 2*E1_brute - S_exact, 5))
print("Sc:", nstr(Sc,25)); print("E1':", nstr(E1_brute,25), " tailM:", nstr(tailM, 12))
'''
exec(code)
