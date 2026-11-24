import mpmath as mp
mp.mp.dps = 60
z = mp.zeta; pi = mp.pi

# (1) classical alternating Euler sum A_s = sum (-1)^{l-1} H_l / l^s, check closed form for even s=2m:
# 2 A_{2m} = (2m+1) eta(2m+1) - zeta(2m+1) - 2 sum_{j=1}^{m-1} eta(2j) zeta(2m+1-2j)
def eta(s): return (1-mp.power(2,1-s))*z(s)
def Asum(s, N=2000000):
    # mpmath summation with extrapolation
    return mp.nsum(lambda l: (-1)**(l-1)*mp.nsum(lambda i: 1/mp.mpf(i), [1, l], method='direct')/l**s,
                   [1, mp.inf], method='richardson+shanks', tol=mp.mpf('1e-55'))
print("A_{2m}: series vs closed form")
for m in range(1, 7):
    s = 2*m
    A = Asum(s)
    cf = ((s+1)*eta(s+1) - z(s+1))/2 - mp.fsum(eta(2*j)*z(s+1-2*j) for j in range(1, m))
    print(f"m={m}: |A_{s} - closed| = {mp.nstr(abs(A-cf),3)}")

# (2) top-generating-function:  U(kap,tau) = sum_n top(Q_n) tau^n
# form A: tau * Re sum_{l>=1} binom(tau,l)/(l - tau(1/2 + i kap))
def U_series(kap, tau):
    return tau*mp.re(mp.nsum(lambda l: mp.binomial(tau, l)/(l - tau*(mp.mpf('0.5')+1j*kap)), [1, mp.inf]))
# form B: even part sinh + odd part digamma
def U_closed(kap, tau):
    u = kap*tau
    sh = pi*u/mp.sinh(pi*u)
    R = mp.re(mp.digamma(1+1j*u)) + mp.euler
    return (1/kap**2)*(1 - sh) + (tau/kap**2)*sh*R

def tops(kap, tau):  # direct from polynomial data
    out = mp.mpf(0)
    for n in range(2, 10):
        r = (n-2)//2
        if n % 2 == 0:
            top = (-1)**(n//2-1)*2*eta(n)
        else:
            m = (n-1)//2
            top = (-1)**(m-1)*(z(n) + 2*mp.fsum(eta(2*j)*z(n-2*j) for j in range(1, m)))
        out += top * kap**(2*r) * tau**n
    return out

print("\nU(kap,tau): series vs closed form vs polynomial tops (n<=9 + O(tau^10))")
for kap, tau in [(mp.mpf(2), mp.mpf('0.04')), (mp.mpf(6), mp.mpf('0.03')), (mp.mpf(1), mp.mpf('0.05'))]:
    a, b, c = U_series(kap,tau), U_closed(kap,tau), tops(kap,tau)
    print(f"kap={kap}, tau={tau}: series={mp.nstr(a,14)}")
    print(f"   closed={mp.nstr(b,14)}  diff_series={mp.nstr(abs(a-b),3)}")
    print(f"   tops  ={mp.nstr(c,14)}  (expected to differ by O(tau^10) ~ subleading top)  diff={mp.nstr(abs(a-c),3)}")
