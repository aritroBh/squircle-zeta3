import mpmath
mpmath.mp.dps = 60
from mpmath import mpf, sqrt, log, quad
from math import comb
S = lambda u: sqrt(1+u**2)
def m_even(u): return (1 - sqrt(2)/S(u))/(1-u**2)
def L1(u): return log((1+u)/2)
def L2(u): return log((1-u)/2)
def Tj(d,j,measure):
    return (-1)**(j+1)*comb(d,j)*quad(lambda u: L1(u)**(d-j)*L2(u)**j*measure(u), [-1,0,1])
for d in [2,4]:
    total = sum(Tj(d,j,m_even) for j in range(0,d//2))
    mid=(d)//2
    total += (-1)**(mid+1)*comb(d,mid)*quad(lambda u: L1(u)**(d-mid)*L2(u)**mid*(1-u)/2*m_even(u), [-1,0,1])
    F_even = -quad(lambda u: log((1+u)/(1-u))**d*m_even(u), [-1,0,1])
    corr_new = -mpf('0.25')*(-1)**(d//2+1)*comb(d,d//2)*quad(lambda u: (L1(u)*L2(u))**(d//2)*m_even(u), [-1,0,1])
    print(f"d={d}: total={mpmath.nstr(total,30)}")
    print(f"   F_even/2={mpmath.nstr(F_even/2,30)}  F_even/2+corr_new={mpmath.nstr(F_even/2+corr_new,30)}")
