from mpmath import mp, mpf, zeta, log, pi, nprint, sin
mp.dps = 50
def psum(f, N, H=False):
    S = mpf(0); Hk = mpf(0)
    for m in range(1, N+1):
        if H: Hk += 1/mpf(m)
        S += f(m, Hk)
    return S, Hk
# S = sum sin(m pi/3) H_m / m^2
for N in [1000, 10000, 100000]:
    S,_ = psum(lambda m,H: sin(m*pi/3)*H/m**2, N, H=True)
    print(f"S (N={N}) =", nprint(S, 25))
print("pi^3/36 =", nprint(pi**3/36, 25))
print("41 pi^3/1296 =", nprint(41*pi**3/1296, 25))
# sum sin(m pi/3)/m^3
for N in [10000, 100000]:
    T,_ = psum(lambda m,H: sin(m*pi/3)/m**3, N)
    print(f"sum sin/ m^3 (N={N}) =", nprint(T, 25))
print("5 pi^3/162 =", nprint(5*pi**3/162, 25))
# sum sin(2m pi/3)/m^3
for N in [100000]:
    T2,_ = psum(lambda m,H: sin(2*m*pi/3)/m**3, N)
    print(f"sum sin(2mpi/3)/m^3 =", nprint(T2, 25))
print("2 pi^3/81 =", nprint(2*pi**3/81, 25))
