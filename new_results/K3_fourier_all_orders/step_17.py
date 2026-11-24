import mpmath as mp
mp.mp.dps = 60
z = mp.zeta; pi = mp.pi

def U_series(kap, tau):
    return 2*tau*mp.re(mp.nsum(lambda l: mp.binomial(tau, l)/(l - tau*(mp.mpf('0.5')+1j*kap)), [1, mp.inf]))
def U_closed(kap, tau):
    u = kap*tau
    sh = pi*u/mp.sinh(pi*u)
    R = mp.re(mp.digamma(1+1j*u)) + mp.euler
    return (1/kap**2)*(1 - sh) + (tau/kap**2)*sh*R
def U_integral(kap, tau):
    return 2*tau*mp.re(mp.quad(lambda s: (mp.exp(tau*mp.log1p(mp.e**(-s)))-1)*mp.exp(tau*s*(mp.mpf('0.5')+1j*kap)), [0, 2, 8, 30, mp.inf]))

print("U: series vs closed vs integral")
for kap, tau in [(mp.mpf(2), mp.mpf('0.04')), (mp.mpf(6), mp.mpf('0.03')), (mp.mpf(1), mp.mpf('0.05')), (mp.mpf(20), mp.mpf('0.02'))]:
    a, b, c = U_series(kap,tau), U_closed(kap,tau), U_integral(kap,tau)
    print(f"kap={kap}, tau={tau}: |series-closed|={mp.nstr(abs(a-b),3)}, |series-integral|={mp.nstr(abs(a-c),3)}")
    print(f"    U = {mp.nstr(a,16)}")

# asymptotic cross-check: U should equal sum over n of top(Q_n) tau^n; here tau small, all n
def tops_all(kap, tau, N=30):
    out = mp.mpf(0)
    for n in range(2, N):
        r = (n-2)//2
        if n % 2 == 0:
            top = (-1)**(n//2-1)*2*(1-mp.power(2,1-n))*z(n)
        else:
            m = (n-1)//2
            top = (-1)**(m-1)*(z(n) + 2*mp.fsum((1-mp.power(2,1-2*j))*z(2*j)*z(n-2*j) for j in range(1, m)))
        out += top * kap**(2*r) * tau**n
    return out
print("\nU_closed vs sum of top(Q_n), n<30, tau=0.01, kap=2,6:")
for kap in [mp.mpf(2), mp.mpf(6)]:
    print(mp.nstr(abs(U_closed(kap, mp.mpf('0.01')) - tops_all(kap, mp.mpf('0.01'))), 3))
