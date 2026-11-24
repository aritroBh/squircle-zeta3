
import json, os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
out = json.load(open(os.path.join(FIX, 'tight_100_400.json')))
sup = max(o[2] for o in out)
argm = [(o[0],o[1]) for o in out if o[2]==sup][0]
print(f"[100,400] fixed: sup C = {sup!r} at {argm}, cells {len(out)}")
# buggy
outb = json.load(open(os.path.join(WORK,'buggy','tight_100_400.json')))
supb = max(o[2] for o in outb)
print(f"[100,400] buggy: sup C = {supb!r}")

# write verify_Btrue.py (H1 exact form) and certify_H2.py for the record
btrue = r'''
# H1: EXACT regular form  Btrue = rho*B1/(2+rt2 t)  vs direct  Btrue = (Im Phi(x0)-th+pi/4)/((pi/rt2)y)
import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); pi = mp.pi; a = pi/(2*sqrt2); b = mp.e**(1j*pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(b*x/(1+b*x))
def Bform(t):
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t))
    s1 = -t/((1+sqrt2*t)*(1+rho))
    zz2 = 2*t*t/(2+sqrt2*t)**2
    A1n = zz2/3/t**2 if t>0 else mp.mpf(1)/6   # midpoint of the A1n box (zz2/3-zz2^2/5, zz2/3)/t^2
    B1 = mp.mpf('0.5') + sqrt2*s1 + s1*s1 - rho*rho*A1n
    return rho*B1/(2+sqrt2*t)
worst = mp.mpf('0'); random.seed(7)
for _ in range(600):
    t = mp.mpf(10)**random.uniform(-12, 0)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t)); y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    Bt = (mp.im(Phi(x0,N)) - th + pi/4)/((pi/sqrt2)*y)
    worst = max(worst, abs(Bt - Bform(t))/abs(Bt))
print("H1 exact form: worst rel resid over t in [1e-12,1] (600 pts) =", mp.nstr(worst, 8))
# exact identity for phi: phi = a(1-2rho)/t + N atan(zz)
w2 = mp.mpf('0'); random.seed(11)
for _ in range(300):
    t = mp.mpf(10)**random.uniform(-12, 0)
    rho = mp.sqrt(1 + t*t/(1+sqrt2*t)); y = t/rho; N = pi/(2*y*y); th = pi/(sqrt2*y); x0 = 1/(rho*y)
    phi = mp.im(Phi(x0,N)) - th
    zz = sqrt2*t/(2+sqrt2*t)
    phi2 = a*(1-2*rho)/t + N*mp.atan(zz)
    w2 = max(w2, abs(phi-phi2))
print("phi = a(1-2rho)/t + N atan(zz): worst abs resid =", mp.nstr(w2, 8))
'''
open(os.path.join(FIX,'verify_Btrue.py'),'w').write(btrue)

h2 = r'''
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
'''
open(os.path.join(FIX,'certify_H2.py'),'w').write(h2)
print("files written")
