
import mpmath as mp
mp.mp.dps = 50
sqrt2 = mp.sqrt(2); pi = mp.pi; sqpi = mp.sqrt(pi)
b = mp.e**(1j*pi/4)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def gfun(x): return Hc(b*x)/x

# Lemma H2a: x^2 |g(x)| <= 1 on (0,1]  (proved via series; verify numerically)
worst = mp.mpf('0'); argx = None
for k in range(1, 20001):
    x = (k/20000)**2   # dense near 0
    v = x*x*abs(gfun(x))
    if v > worst: worst, argx = v, x
print("H2a numeric: max x^2|g(x)| over grid =", mp.nstr(worst, 15), "at x =", argx, "(<=1, limit at 0 is 1)")

# Series check: sHc(s) = 1 - sum c_k s/(s+2k), c_k = (-1)^{k+1} binom(1/2,k) > 0, sum c_k = 1
from sympy import binomial, Rational, summation, symbols, simplify
k = symbols('k', positive=True, integer=True)
csum = summation((-1)**(k+1)*binomial(Rational(1,2),k), (k,1,10))
print("partial sum c_k (k=1..10) =", csum, "-> 1 as K->inf since (1-u)^{1/2}->0 at u=1")

# numeric check of series identity sHc(s) = 1 - sum_{k>=1} c_k s/(s+2k) at a few s on the ray
def c_k(k): return (-1)**(k+1)*mp.binomial(mp.mpf('0.5'), k)
for xv in ['0.13','0.5','1.0','3.7']:
    s = mp.mpf(xv)*b
    lhs = s*Hc(s)
    rhs = 1 - mp.nsum(lambda kk: c_k(int(kk))*s/(s+2*kk), [1, mp.inf])
    print(f"s={xv}*b: |sHc(s)|={mp.nstr(abs(lhs),12)}  series resid={mp.nstr(abs(lhs-rhs),6)}  |rhs|={mp.nstr(abs(rhs),12)}")

# key inequality used in proof: 2 Re w >= |w| for w = s/(s+2k), s = x e^{i pi/4}
worst2 = mp.mpf('9')
for kx in range(1, 50):
    for xv in [mp.mpf(10)**(mp.mpf(j)/10) for j in range(-30, 21)]:
        s = xv*b; w = s/(s+2*kx)
        r = 2*mp.re(w)/abs(w)
        worst2 = min(worst2, r)
print("min 2Re(w)/|w| over grid =", mp.nstr(worst2, 10), " (>1 => 2Re w_k >= |w_k| confirmed)")

# Lemma H2b: 1+sqrt2 x+x^2 >= (2+sqrt2) x on [0,1]
diff = lambda x: 1+sqrt2*x+x**2 - (2+sqrt2)*x
print("H2b: (1-x)(x-... ) check: min over grid =", min(diff(k/1000) for k in range(1001)), " = x+1/x-2 form >=0 by AM-GM")

# TL0 consequence: I0 = int_0^1 e^{RePhi}|g| dx <= 2 (2+sqrt2)^{-N/2}/(N-2), and the e^theta-unscaled TL0
a = pi/(2*sqrt2)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))
for Nv in [100, 400, 1000, 10000]:
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    f = lambda x: abs(mp.e**Phi(x,N)*gfun(x))
    I0 = mp.quad(f, [0, mp.mpf('0.5'), mp.mpf('0.95'), 1])
    proved = 2*(2+sqrt2)**(-N/2)/(N-2)
    coded = sqrt2*(2+sqrt2)**(-N/2)/(N-1)  # coded TL0 with e^theta removed
    print(f"N={Nv}: I0 = {mp.nstr(I0,8)}  proved_bd = {mp.nstr(proved,8)}  ratio I0/proved = {mp.nstr(I0/proved,6)}  I0/coded = {mp.nstr(I0/coded,6)}")
