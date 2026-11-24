import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

# sup x^2 |g(x)| on [0,1]
smax = mp.mpf('0'); argm=None
for k in range(0, 2001):
    x = (k/2000)**2   # denser near 0
    v = x*x*abs(g(x)) if x>0 else mp.mpf('0')
    if v > smax: smax, argm = v, x
print("sup x^2|g(x)| on [0,1]:", mp.nstr(smax,12), "near x =", argm, " (limit x->0: 1)")
print("sqrt2/2 =", mp.nstr(sqrt2/2,10), " sqrt2 =", mp.nstr(sqrt2,10))

for Nv in [100, 400, 1000]:
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    cc0 = (2+sqrt2)/4
    # I0 = ∫_0^1 |e^{Phi(x)} g(x)| dx ; integrand ~ x^{N-2} well-behaved for large N
    f = lambda x: abs(mp.e**Phi(x,N)*g(x))
    # split [0, 1/2] and [1/2, 1] for safety
    I0 = mp.quad(f, [0, mp.mpf('0.5'), mp.mpf('0.9'), 1])
    ex0 = mp.log(1/(2*mp.sqrt(cc0))) + mp.sqrt(pi/N)
    TL0 = sqrt2*mp.e**(ex0*N)/(N-1)
    print(f"N={Nv}: I0*e^th = {mp.nstr(I0*mp.e**th,10)}  TL0 = {mp.nstr(TL0,10)}  ratio I0 e^th/TL0 = {mp.nstr(I0*mp.e**th/TL0,10)}")
    # my derived bound: 2*supx2g*(2 sqrt(cc0))^{-N}/(N-1)
    mine = 2*smax*(2*mp.sqrt(cc0))**(-N)/(N-1)
    print(f"     my bound = {mp.nstr(mine,10)}  (mine/TL0 = {mp.nstr(mine/TL0,6)})")
