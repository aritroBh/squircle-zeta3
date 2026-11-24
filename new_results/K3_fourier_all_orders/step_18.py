import mpmath as mp
mp.mp.dps = 50
pi = mp.pi; s2 = mp.sqrt(2); z = mp.zeta

def h(tau, t):
    q = 1/tau
    return (abs(mp.cos(t))**q + abs(mp.sin(t))**q)**tau
def h1(t): return max(abs(mp.cos(t)), abs(mp.sin(t)))

# (a) exact boundary-layer identity: h_p - h_1 = (1+e^{-2y})^{-1/2} (e^{tau L(y/tau)} - 1), y = -log tan theta on [0,pi/4]
print("(a) exact profile identity check")
for tau in [mp.mpf('0.1'), mp.mpf('0.02')]:
    for th in [pi/16, pi/8, pi/4 - 3*tau, pi/4 - mp.mpf('0.5')*tau, pi/4]:
        y = -mp.log(mp.tan(th))
        prof = (1+mp.e**(-2*y))**mp.mpf('-0.5')*(mp.exp(tau*mp.log1p(mp.e**(-y/tau)))-1)
        print(f"  tau={tau}, th={mp.nstr(th,6)}: diff = {mp.nstr(abs((h(tau,th)-h1(th))-prof),3)}")

# (b) partial Fourier sums reproduce h_p - h_1; and profile (tau/sqrt2) log(1+e^{-s}) in the layer
def a_minus_A(k, tau):
    Gk = lambda x: (lambda th_: mp.cos(4*k*th_)*mp.cos(th_)**2*mp.sin(th_))(mp.atan(mp.e**(-x)))
    return 8/pi*tau*mp.quad(lambda s: (mp.e**(tau*mp.log1p(mp.e**(-s)))-1)*Gk(tau*s), [0, 2, 8, 30, mp.inf])
def A(k): return -(-1)**k*4*s2/(pi*(16*k*k-1))

tau = mp.mpf('0.02')
K = 120
aks = [A(k)+a_minus_A(k, tau) for k in range(K+1)]
def fourier_sum(th):
    return aks[0]/2 + mp.fsum(aks[k]*mp.cos(4*k*th) for k in range(1, K+1))
print("\n(b) partial Fourier sum (K=120) vs h_p-h_1, tau=0.02")
for th in [pi/16, pi/8, pi/4-0.1, pi/4-0.01, pi/4-0.002, pi/4]:
    err = fourier_sum(th) - (h(tau,th)-h1(th))
    print(f"  th={mp.nstr(th,8)}: sum-h = {mp.nstr(err,3)}")
print("\n  normalized: (sum - (h_p-h_1)) small; (h_p-h_1)/(tau/sqrt2) vs log(1+e^{-s}), s=-log(tan th)/tau")
for th in [pi/4-0.08, pi/4-0.04, pi/4-0.02, pi/4-0.005]:
    s = -mp.log(mp.tan(th))/tau
    print(f"  th=pi/4-{mp.nstr(pi/4-th,6)}: s={mp.nstr(s,6)}, Fourier/(tau/sqrt2)={mp.nstr(fourier_sum(th)/(tau/s2),10)}, log(1+e^-s)={mp.nstr(mp.log1p(mp.e**(-s)),10)}")
