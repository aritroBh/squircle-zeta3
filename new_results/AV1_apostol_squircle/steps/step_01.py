
# step_01.py : A(p) = iint_{D_p} dx dy/(1-xy), two independent methods, 110 digits.
# Method A: 1D integral  int_0^1 -log(1 - x*(1-x^p)^(1/p))/x dx   (mpmath quad)
# Method B: Dirichlet series  sum_{m>=1} Gamma(1+m/p)^2 / (m^2 Gamma(1+2m/p))
#           (from 1/(1-xy) = sum (xy)^n and Dirichlet's integral on the squircle)
import mpmath as mp
mp.mp.dps = 110

def A_quad(p):
    p = mp.mpf(p)
    def yp(x):
        if x == 0: return mp.mpf(1)
        return (1 - x**p)**(1/p)
    f = lambda x: -mp.log(1 - x*yp(x))/x if x != 0 else mp.mpf(1)
    # split at points where integrand's derivatives are tricky: near x=1 drop like (1-x)^(1/p)
    return mp.quad(f, [0, mp.mpf('0.5'), mp.mpf('0.9'), mp.mpf('0.99'), 1-1e-6, 1])

def A_series(p, tol=mp.mpf('1e-105')):
    p = mp.mpf(p)
    s = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma(1 + m/p)**2 / (m**2 * mp.gamma(1 + 2*m/p))
        s += t
        if m > 20 and abs(t) < tol: break
    return s

# sanity: Dirichlet formula area(D_p) = Gamma(1+1/p)^2/Gamma(1+2/p); at p=2 must be pi/4
assert abs(mp.gamma(1+mp.mpf('0.5'))**2/mp.gamma(2) - mp.pi/4) < mp.mpf('1e-100')

ps = ['0.5','0.75','1','1.25','1.5','2','3','4','8','16','64','256']
print(f"{'p':>6} | {'A_quad':>30} | {'A_series':>30} | {'|diff|':>8}")
for p in ps:
    a = A_quad(p); b = A_series(p)
    d = abs(a-b)
    print(f"{p:>6} | {mp.nstr(a,25):>30} | {mp.nstr(b,25):>30} | {mp.nstr(d,3):>8}")

print()
print('zeta(2) =', mp.nstr(mp.zeta(2), 30))
print('A(256) - zeta(2) =', mp.nstr(A_series(256) - mp.zeta(2), 10))
