
# H2: the two lemmas behind TL0 = 2 (2+rt2)^{-N/2} e^{theta}/(N-2)
import mpmath as mp
mp.mp.dps = 50
sqrt2 = mp.sqrt(2); pi = mp.pi; sqpi = mp.sqrt(pi); b = mp.e**(1j*pi/4)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(b*x)/x
# (i) x^2|g(x)| <= 1 on (0,1]  [proof: sHc(s) = 1 - sum c_k w_k, |.|^2 = 1 - 2A + |sum c w|^2 <= 1 since 2Re w_k >= |w_k|]
mx = max((k/20000)**4*abs(g((k/20000)**2)) for k in range(1, 20001))
print("(i) max x^2|g(x)| grid =", mp.nstr(mx, 12), "(<=1 confirmed; sup = 1 as x->0)")
# series identity check
def c_k(k): return (-1)**(k+1)*mp.binomial(mp.mpf('0.5'), k)
s = mp.mpf('0.37')*b
res = s*Hc(s) - (1 - mp.nsum(lambda kk: c_k(int(kk))*s/(s+2*kk), [1, mp.inf]))
print("    series identity resid at s=0.37b:", mp.nstr(abs(res), 8))
wmin = min(2*mp.re((b*mp.mpf(k)/10)/(b*mp.mpf(k)/10+2*j))/abs((b*mp.mpf(k)/10)/(b*mp.mpf(k)/10+2*j))
           for k in range(1, 200) for j in range(1, 40))
print("    min 2Re(w)/|w| =", mp.nstr(wmin, 10), "(>=1 confirms 2Re w_k >= |w_k|; limit sqrt2 as k->inf)")
# (ii) AM-GM: 1+sqrt2 x+x^2 >= (2+sqrt2) x on [0,1]
print("(ii) min (1+rt2 x+x^2) - (2+rt2)x on [0,1] =", min(1+sqrt2*k/1000+(k/1000)**2-(2+sqrt2)*k/1000 for k in range(1001)))
# (iii) the bound vs true I0
a = pi/(2*sqrt2)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))
for Nv in [100, 400, 1000]:
    N = mp.mpf(Nv)
    I0 = mp.quad(lambda x: abs(mp.e**Phi(x,N)*g(x)), [0, mp.mpf('0.5'), mp.mpf('0.95'), 1])
    bd = 2*(2+sqrt2)**(-N/2)/(N-2)
    print(f"    N={Nv}: I0 = {mp.nstr(I0,8)} <= 2(2+rt2)^{{-N/2}}/(N-2) = {mp.nstr(bd,8)}  ratio {mp.nstr(I0/bd,6)}")
