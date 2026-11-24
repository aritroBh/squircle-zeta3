import mpmath as mp, time
mp.mp.dps = 450

def Li21_fast(a, c, digits=450):
    a, c = mp.mpf(a), mp.mpf(c)
    # terms a^n c^m /(n^2 m), n>m>=1; tail estimate
    N = 2
    while True:
        # worst term at n=N+1,m=N: |a|^(N+1)|c|^N
        if abs(a)**(N+1)*abs(c)**N < mp.mpf(10)**(-digits-20) and abs(a)<=1:
            break
        N += 1
    apw = [mp.mpf(1)]; cpw = [mp.mpf(1)]
    for i in range(N+2):
        apw.append(apw[-1]*a); cpw.append(cpw[-1]*c)
    tot = mp.mpf(0)
    for n in range(2, N+2):
        an = apw[n]; inv = 1/(mp.mpf(n)*mp.mpf(n))
        s = mp.fsum([an*cpw[m]/m for m in range(1,n)])*inv
        tot += s
    return tot, N

t0=time.time()
x1 = mp.mpf(-1)/3; x2 = mp.mpf(-1)/9
Lam1, N1 = Li21_fast(x1, x2/x1)
t1=time.time()
print(f"Lam1 (N={N1}): {mp.nstr(Lam1, 30)}  [{t1-t0:.1f}s]")
# benchmark pslq at dim ~60, 450 dps
import random
basis = [mp.mpf(1), mp.log(2), mp.log(3), mp.log(5), mp.pi**2/6, mp.zeta(3), Lam1]
for i in range(50):
    basis.append(mp.log(2)**random.randint(0,3)*mp.log(3)**random.randint(0,2)*mp.mpf(random.randint(1,5)))
t2=time.time()
rel = mp.pslq(basis, tol=mp.mpf('1e-440'), maxcoeff=10**20, maxsteps=5000)
t3=time.time()
print(f"pslq dim {len(basis)} at 450dps: {t3-t2:.1f}s, result: {rel}")
