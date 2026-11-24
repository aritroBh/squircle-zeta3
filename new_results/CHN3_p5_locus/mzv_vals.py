"""Numerical double zeta values zeta(a,b) = sum_{n>m>=1} n^-a m^-b  (a>=2, b>=1).
Method: zeta(a,b) = sum_{n>=1} n^-a H^{(b)}_{n-1}. Sum n < N directly; for the tail use
H^{(b)}_{n-1} = zeta(b) - zeta(b,n) (b>=2) and the Euler-Maclaurin expansion of the Hurwitz zeta
zeta(b,n) ~ n^{1-b}/(b-1) + n^{-b}/2 + sum_j B_{2j}/(2j)! (b)_{2j-1} n^{-b-2j+1}, so that each tail
piece is a Hurwitz zeta zeta(s, N). (b=1 handled with digamma expansion.)"""
from mpmath import mp, mpf, zeta, bernoulli, rf, factorial, harmonic, euler, log
def dz(a, b, N=400, J=40):
    a = int(a); b = int(b)
    S = mpf(0); H = mpf(0)
    for n in range(1, N):
        S += H / mpf(n)**a
        H += mpf(1) / mpf(n)**b
    if b >= 2:
        tail = zeta(b) * zeta(a, N) - zeta(a + b - 1, N) / (b - 1) - zeta(a + b, N) / 2
        for j in range(1, J):
            tail -= bernoulli(2*j) / factorial(2*j) * rf(b, 2*j - 1) * zeta(a + b + 2*j - 1, N)
    else:  # H_{n-1} = log n + gamma - 1/(2n) - sum_j B_{2j}/(2j) n^{-2j}
        from mpmath import diff
        tail = euler * zeta(a, N) - zeta(a, N, 1) - zeta(a + 1, N) / 2  # -d/ds zeta(s,N) = sum log n / n^s
        for j in range(1, J):
            tail -= bernoulli(2*j) / (2*j) * zeta(a + 2*j, N)
    return S + tail
if __name__ == "__main__":
    mp.dps = 110
    z53 = dz(5, 3); z35 = dz(3, 5)
    print("zeta(5,3) =", mp.nstr(z53, 90))
    print("zeta(3,5) =", mp.nstr(z35, 90))
    print("stuffle check z53+z35+z8-z3z5:", mp.nstr(z53 + z35 + zeta(8) - zeta(3)*zeta(5), 5))
    print("zeta(2,1)-zeta(3):", mp.nstr(dz(2, 1) - zeta(3), 5))
    print("zeta(3,1)-zeta(4)/4:", mp.nstr(dz(3, 1) - zeta(4)/4, 5))
    print("zeta(4,2)-(zeta(3)^2-4/3 zeta(6)):", mp.nstr(dz(4, 2) - (zeta(3)**2 - mpf(4)/3*zeta(6)), 5))
