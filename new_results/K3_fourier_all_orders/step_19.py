import mpmath as mp
mp.mp.dps = 50
pi = mp.pi; s2 = mp.sqrt(2)

def h(tau, t):
    q = 1/tau
    return (abs(mp.cos(t))**q + abs(mp.sin(t))**q)**tau
def h1(t): return max(abs(mp.cos(t)), abs(mp.sin(t)))
def a_minus_A(k, tau):
    Gk = lambda x: (lambda th_: mp.cos(4*k*th_)*mp.cos(th_)**2*mp.sin(th_))(mp.atan(mp.e**(-x)))
    return 8/pi*tau*mp.quad(lambda s: (mp.e**(tau*mp.log1p(mp.e**(-s)))-1)*Gk(tau*s), [0, 2, 8, 30, mp.inf])
def A(k): return -(-1)**k*4*s2/(pi*(16*k*k-1))

tau = mp.mpf('0.02')
K = 200
dk = [a_minus_A(k, tau) for k in range(K+1)]  # a_k - A_k
def diff_sum(th):
    return dk[0]/2 + mp.fsum(dk[k]*mp.cos(4*k*th) for k in range(1, K+1))
print("(b) difference-series (K=200) vs h_p - h_1, tau=0.02")
for th in [pi/16, pi/8, pi/4-0.1, pi/4-0.01, pi/4]:
    print(f"  th={mp.nstr(th,8)}: err = {mp.nstr(diff_sum(th)-(h(tau,th)-h1(th)),3)}  (h-h1={mp.nstr(h(tau,th)-h1(th),8)})")
print("\n  boundary layer: (sum)/(tau/sqrt2) vs log(1+e^{-s})")
for th in [pi/4-0.08, pi/4-0.04, pi/4-0.02, pi/4-0.005, pi/4]:
    s = -mp.log(mp.tan(th))/tau
    print(f"  th=pi/4-{mp.nstr(pi/4-th,8)}: s={mp.nstr(s,8)}, sum/(tau/sqrt2)={mp.nstr(diff_sum(th)/(tau/s2),12)}, log(1+e^-s)={mp.nstr(mp.log1p(mp.e**(-s)),12)}")
