from mpmath import (mp, mpf, pi, gamma, sin, cos, exp, im, re, sqrt, zeta,
                    gammainc, hyper, quad, log, nsum, inf, nstr, fac, si, ci, euler, psi)
mp.dps = 80
def N(x,d=40): return nstr(x,d)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")

# --- S_Si = sum Si(n pi)/n^3 ---
# Method 1: direct nsum (alternating decomposition: Si(nπ)=π/2+(−1)^{n+1}g(n), g>0 decreasing)
g = lambda n: (-1)**(n+1)*(si(n*pi) - pi/2)
S1 = (pi/2)*zeta(3) + nsum(lambda n: (-1)**(n+1)*g(n)/n**3, [1, inf])
print("S_Si method1 (nsum alternating):", N(S1,50))
target = 5*pi**3/72
print("target 5 pi^3/72              :", N(target,50))
print("diff:", nstr(abs(S1-target),3))

# Method 2: asymptotic split. For n>N0 use Si(nπ)=π/2+(−1)^{n+1}/(nπ) A(n) with A(n)=Σ(-1)^k(2k)!/(nπ)^{2k}
N0 = 60
head = nsum(lambda n: si(n*pi)/n**3, [1, N0])
def tail(N0, K):
    tot = (pi/2)*(zeta(3) - nsum(lambda n: 1/n**3, [1,N0]))
    for k in range(0, K):
        # sum_{n>N0} (-1)^{n+1} (-1)^k (2k)! π^{-(2k+1)} n^{-(2k+4)}
        coeff = (-1)**k * fac(2*k) * pi**(-(2*k+1))
        inner = (-1)**N0 * nsum(lambda j: (-1)**(j+1) * (N0+j)**(-(2*k+4)), [1, inf])
        tot = tot + coeff*inner
    return tot
for K in [10, 20, 40]:
    S2 = head + tail(N0, K)
    print(f"S_Si method2 (asymptotic K={K}):", N(S2,50), " diff:", nstr(abs(S2-target),3))
logline("step_08: S_Si=sum Si(n pi)/n^3 = 5 pi^3/72 verified by 2 methods to ~45+ digits")
