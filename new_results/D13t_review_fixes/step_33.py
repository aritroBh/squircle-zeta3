import mpmath as mp
mp.mp.dps = 50
mpr2=mp.sqrt(2); mpp=mp.pi; sqpi=mp.sqrt(mpp); mbb=mp.e**(1j*mpp/4); ma=mpp/(2*mpr2)
# 2 Re w_k >= |w_k| dense scan, w_k = s/(s+2k), s = b x, x in (0,1], k >= 1
mn = mp.mpf('9'); at=None
xs = [mp.mpf(i)/1000 for i in range(1,1001)] + [mp.mpf(10)**(-j) for j in range(1,7)]
for x in xs:
    s = x*mbb
    for kk in list(range(1,300)) + [10**3, 10**4, 10**5, 10**6, 10**8]:
        w = s/(s+2*kk)
        r = 2*mp.re(w)/abs(w)
        if r < mn: mn=r; at=(x,kk)
print(f"R4: min 2Re(w_k)/|w_k| over dense (x,k) = {mp.nstr(mn,10)} at x={at[0]}, k={at[1]} (>=1 required; theory min -> 1 as x->inf, -> rt2 as k->inf or x->0)")
# I0 vs bound at several N
def Phi(x,N): return (-1+1j)*ma*x + N*mp.log(mbb*x/(1+mbb*x))
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(mbb*x)/x
for Nv in [100, 400, 1000, 10000]:
    N = mp.mpf(Nv)
    I0 = mp.quad(lambda x: abs(mp.e**Phi(x,N)*g(x)), [0, mp.mpf('0.5'), mp.mpf('0.99'), 1])
    bd = 2*(2+mpr2)**(-N/2)/(N-2)
    print(f"  N={Nv}: I0 = {mp.nstr(I0,6)} <= 2(2+rt2)^(-N/2)/(N-2) = {mp.nstr(bd,6)} ; ratio = {mp.nstr(I0/bd,6)}")
# e^{Re Phi(x)} = e^{-a x} x^N (1+rt2 x+x^2)^{-N/2} exactness
N = mp.mpf(777)
wr = max(abs(mp.e**mp.re(Phi(x,N)) - mp.e**(-ma*x)*x**N*(1+mpr2*x+x*x)**(-N/2))/mp.e**mp.re(Phi(x,N)) for x in [mp.mpf(i)/40 for i in range(1,41)])
print(f"  e^RePhi factorisation: worst rel resid on (0,1] = {mp.nstr(wr,6)}")
