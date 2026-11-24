import mpmath as mp
mp.mp.dps = 80
mpr2=mp.sqrt(2); mpp=mp.pi; sqpi=mp.sqrt(mpp); mbb=mp.e**(1j*mpp/4)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def c_k(k): return (-1)**(k+1)*mp.binomial(mp.mpf('0.5'), k)
# sum c_k = 1 via nsum to inf
S = mp.nsum(c_k, [1, mp.inf], tol=mp.mpf('1e-60'))
print("sum_{k>=1} c_k =", mp.nstr(S, 20), " (resid vs 1:", mp.nstr(abs(S-1),6), ")")
# closed form: 1-(1-u)^{1/2} = sum c_k u^k; evaluate at u=0.9876543 as an identity check
u = mp.mpf('0.9876543')
lhs = 1-(1-u)**mp.mpf('0.5'); rhs = mp.nsum(lambda k: c_k(int(k))*u**int(k), [1, mp.inf], tol=mp.mpf('1e-70'))
print("generating-function identity resid at u=0.9876543:", mp.nstr(abs(lhs-rhs),6))
# series identity sHc(s) with controlled truncation: sum to K=2e6 exactly + tail bound
def sHc_series(s, K=200000):
    return 1 - mp.fsum(c_k(k)*s/(s+2*k) for k in range(1,K))
for xr in [mp.mpf('0.001'), mp.mpf('0.37'), mp.mpf('1.0')]:
    s = xr*mbb
    v1 = s*Hc(s); v2 = sHc_series(s)
    print(f"  x={xr}: |sHc(s) - (1-sum_{{k<=2e5}})| = {mp.nstr(abs(v1-v2),6)}")
# beta integral with singularity handled: substitute u=v^2
s = mp.mpf('0.37')*mbb
be = mp.quad(lambda v: 2*v*(v**2)**(s/2-1)*(1-v**2)**mp.mpf('0.5'), [0,1])/2
print("beta-integral identity resid (u=v^2 trick):", mp.nstr(abs(Hc(s)-be),6))
