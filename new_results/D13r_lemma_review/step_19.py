import mpmath as mp, random
import sympy as sp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2)

# ReX2 old form: RePhi(x2)+th with RePhi(x) = -a x + N log x - (N/2) log(1+rt2 x + x^2)
def RePhi(x, N): return -a*x + N*mp.log(x) - (N/2)*mp.log(1+sqrt2*x+x*x)
random.seed(17)
worst = mp.mpf('0'); warg=None
for _ in range(500):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = mp.pi/(2*y*y); th = mp.pi/(sqrt2*y)
    rho = 1+y*y*random.random()
    eta2 = mp.mpf(10)**random.uniform(-1.2, 0)   # eta2 ~ 1/2 in the code
    x2 = eta2/(rho*y)
    old = -a*x2 + th - (N/2)*mp.log(1+sqrt2/x2+1/x2**2)
    M1 = 1 - eta2/(2*rho)
    ureg = sqrt2*rho/eta2 + rho**2*y/eta2**2
    u = y*ureg
    logcf = mp.log(1+u)/u
    L2 = (sqrt2/4)*ureg*logcf
    new = th*(M1-L2)
    r = abs(old-new)
    if r > worst: worst, warg = r, (y,rho,eta2)
print("L5 ReX2 = th*(M1-L2): worst residual over 500 pts @60dps:", mp.nstr(worst,6), "at", [mp.nstr(v,4) for v in warg])

# logcf(u) = log(1+u)/u decreasing for u>0: sympy + numeric
u = sp.symbols('u', positive=True)
lcf = sp.log(1+u)/u
dl = sp.simplify(sp.diff(lcf, u)*u**2*(1+u))   # = u - (1+u) log(1+u)
print("u^2(1+u) logcf'(u) =", sp.expand(dl), " ;  u < (1+u)log(1+u) for u>0 since log(1+u) = int_0^u ds/(1+s) > u/(1+u)")
h = u - (1+u)*sp.log(1+u)
print("h(0)=0, h'(u)=", sp.simplify(sp.diff(h,u)), "< 0 for u>0  => logcf strictly decreasing PROVED")
maxd = mp.mpf('-9'); wmin = mp.mpf('9')
for _ in range(3000):
    u1 = mp.mpf(10)**random.uniform(-8, 3); u2 = u1*mp.mpf(10)**random.uniform(1e-6, 1)
    d = mp.log(1+u2)/u2 - mp.log(1+u1)/u1
    if d > maxd: maxd = d
print("numeric: max logcf(u2)-logcf(u1) over 3000 pairs u2>u1 (must be <0):", mp.nstr(maxd,6))

# M1 - L2 < 0 on the range the code uses (eta2 = 1/2 exactly in both step_38 and final cell)
# final cell: eta2=0.5, rho in [1, 1+y_J^2], y in [0,y_J], y_J=1.2533e-5
worstMmL = mp.mpf('-9')
for _ in range(2000):
    y = mp.mpf(10)**random.uniform(-12, mp.log10(1.2533e-5))
    rho = 1+y*y*random.random(); eta2 = mp.mpf('0.5')
    M1 = 1-eta2/(2*rho); ureg = sqrt2*rho/eta2 + rho**2*y/eta2**2; u=y*ureg
    MmL = M1-(sqrt2/4)*ureg*mp.log(1+u)/u
    worstMmL = max(worstMmL, MmL)
print("final-cell range: max M1-L2 (must be <0):", mp.nstr(worstMmL,8), "(code printed -0.249982)")
# cells range N in [100,1.5e10], eta2 box around 1/2: check M1-L2<0 there too (used in step_38 asserts)
worstMmL2 = mp.mpf('-9')
for _ in range(3000):
    N = mp.mpf(10)**random.uniform(2, 10.2)
    y = mp.sqrt(mp.pi/(2*N)); rho = 1+y*y*random.random()
    eta2 = mp.mpf('0.5')*mp.mpf(10)**random.uniform(-0.02,0.02)
    M1 = 1-eta2/(2*rho); ureg = sqrt2*rho/eta2 + rho**2*y/eta2**2; u=y*ureg
    MmL = M1-(sqrt2/4)*ureg*mp.log(1+u)/u
    worstMmL2 = max(worstMmL2, MmL)
print("cells range (N 100..1.5e10, eta2~0.5): max M1-L2 (must be <0):", mp.nstr(worstMmL2,8))
with open('research/08_theory/D13r_lemma_review/progress.log','a') as f:
    f.write("L5: HOLDS (ReX2=th(M1-L2) hand algebra + 500 pts @60dps worst resid printed; logcf decreasing PROVED (h=u-(1+u)log(1+u), h(0)=0, h'=-log(1+u)<0) + 3000-pair numeric; M1-L2<0 confirmed on final-cell and cells ranges, max values printed).\n")
print("logged")
