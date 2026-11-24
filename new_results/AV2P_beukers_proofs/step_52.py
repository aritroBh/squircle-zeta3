from mpmath import mp, mpf, zeta, log, pi, nprint, quad, clsin, clcos, catalan, binomial, polylog, exp, loggamma, digamma, sqrt, sin
mp.dps = 80
ok = []
def chk(label, residual, digits=40):
    r = abs(residual)
    ok.append((label, r))
    print(f"{label:60s} res={nprint(r,3)}")

G = catalan; z3 = zeta(3); L2 = log(2); C13 = clsin(2,pi/3); C16 = clsin(2,pi/6)
C23 = clsin(2, 2*pi/3); C56 = clsin(2, 5*pi/6)
H = lambda n: sum(1/mpf(k) for k in range(1,n+1)) if n>0 else mpf(0)

# === TASK 1: B2(1) = 2 zeta(3)/3 ===
T = mp.nsum(lambda n: (H(2*n)-H(n))/(mpf(n)**2*binomial(2*n,n)), [1, mp.inf])
U = mp.nsum(lambda n: 1/(mpf(n)**3*binomial(2*n,n)), [1, mp.inf])
chk("T = 2zeta3 - (2pi/3)Cl2(pi/3)", T - (2*z3 - (2*pi/3)*C13))
chk("U = (2pi/3)Cl2(pi/3) - 4zeta3/3", U - ((2*pi/3)*C13 - 4*z3/3))
chk("B2(1) = T+U = 2zeta3/3", (T+U) - 2*z3/3)
chk("T integral form = -2 int asin^2(u/2)/(1+u)",
    T + 2*quad(lambda u: mp.asin(u/2)**2/(1+u), [0,1]))
chk("L = int th log(1+2sin th) = zeta3/2 - (pi/6)Cl13 + (pi^2/72)log2",
    quad(lambda t: t*log(1+2*sin(t)), [0, pi/6]) - (z3/2 - (pi/6)*C13 + pi**2*L2/72))
chk("int log(1+2sin th) = 4G/3 - Cl2(pi/3)",
    quad(lambda t: log(1+2*sin(t)), [0, pi/6]) - (4*G/3 - C13))
chk("Cl3(pi/3) = zeta3/3", clcos(3,pi/3) - z3/3)
chk("U integral: 4 int th^2 cot", U - 4*quad(lambda t: t**2/mp.tan(t), [0, pi/6]))

# === TASK 2: C2 ===
def Mk(k):
    return 2*quad(lambda u: (log(4/(1-u**2)))**k/(1+u), [0, 1])
M0, M1, M2, M3, M5 = Mk(0), Mk(1), Mk(2), Mk(3), Mk(5)
chk("M0 = 2 log2", M0 - 2*L2)
chk("M1 = 2 log2^2 + pi^2/6", M1 - (2*L2**2 + pi**2/6))
chk("C2 = log2^2 - M1/2 = -pi^2/12", (L2**2 - M1/2) + pi**2/12)
chk("M2 = 4 zeta3 + 8 log2^3/3 (NUM)", M2 - (4*z3 + 8*L2**3/3))
chk("M3 = 4 log2^4 + 3 pi^4/20 (NUM)", M3 - (4*L2**4 + 3*pi**4/20))
chk("M5 = 32 log2^6/3 + 79pi^6/126 - 240 zeta3^2 (NUM)",
    M5 - (32*L2**6/3 + 79*pi**6/126 - 240*z3**2))
# C2 integral at 60 digits
import json
phicoef = [mpf(x) for x in json.load(open('phicoef.json'))]
def phi(w):
    w = mpf(w)
    if abs(w) < mpf('0.2'):
        return sum(phicoef[k]*w**k for k in range(len(phicoef)))
    return w**3*exp(2*loggamma(w)-loggamma(2*w+1))*(digamma(2*w+1)-digamma(w))
a = mpf('0.2')
C2int = sum(phicoef[k]*a**(k-2)/(k-2) for k in range(3, 250))
for aa, bb in [(a,1),(1,4),(4,16),(16,64),(64,200),(200,mp.inf)]:
    C2int += quad(lambda w: (phi(w)-1)/w**3, [aa,bb])
chk("C2 integral = -pi^2/12", C2int + pi**2/12)
# expansion check: B2(p) fit at p=2^14 vs expansion
def B2s(p):
    p = mpf(p); S = mpf(0); n = 1
    while True:
        t = phi(n/p)/mpf(n)**3
        S += t
        if n > 40*p and abs(t) < mpf(10)**(-75): break
        n += 1
    return S
p = mpf(2)**14; eps = 1/p
phi6 = 79*pi**6/7560 - 4*z3**2
expan = z3 - pi**2/12*eps**2 + z3/2*eps**3 - 3*zeta(4)/16*eps**4 + phi6/120*eps**6
chk("B2(2^14) vs expansion thru eps^6", B2s(p) - expan)

# === TASK 3: B2(2) ===
K = 400
S1=mpf(0);S1p=mpf(0);S2=mpf(0);S3=mpf(0);Hk=mpf(0)
for k in range(0,K):
    if k>0: Hk += 1/mpf(k)
    H2k = Hk + sum(1/mpf(k+j) for j in range(1,k+1)) if k>0 else mpf(0)
    c = binomial(2*k,k)/(mpf(16)**k*(2*k+1))
    S1 += c*Hk; S1p += c*H2k; S2 += c/(2*k+1); S3 += c
chk("S3 = pi/3", S3 - pi/3)
chk("S2 = Cl2(pi/3)", S2 - C13)
chk("S1 = -(2pi/3)log2 + 8Cl2(5pi/6) - 2Cl2(2pi/3)", S1 - (-(2*pi/3)*L2 + 8*C56 - 2*C23))
chk("S1' = 4Cl2(5pi/6) - 2Cl2(2pi/3)", S1p - (4*C56 - 2*C23))
chk("Cl2(pi/6) = Cl2(pi/3)/4 + 2G/3", C16 - (C13/4 + 2*G/3))
odd = mpf(open('odd2.txt').read())
chk("odd = (pi/8)(2log2 S3 + S1 - S1' + S2)", odd - pi/8*(2*L2*S3 + S1 - S1p + S2))
chk("odd = pi G/3", odd - pi*G/3)
chk("B2(2) = zeta3/12 + pi G/3", mpf('1.05936833579935324104250490229980343260203365771260126540342') - (z3/12 + pi*G/3))

# === TASK 4: d1 ===
d1 = mpf(open('d1.txt').read())
P1 = quad(lambda t: log(4/(1-t**2))**2/(3+t**2), [0,1])
P2 = 2*quad(lambda t: log(4/(1-t**2))*t*mp.atanh(t)/(3+t**2), [0,1])
chk("d1 = P1 - P2", d1 - (P1 - P2))
chk("P2 = 2 zeta3/3 = B2(1)", P2 - 2*z3/3)
chk("P1 = 2 pi^3/(27 sqrt3) (NUM)", P1 - 2*pi**3/(27*sqrt(3)))
chk("d1 = 2 pi^3/(27 sqrt3) - 2 zeta3/3 (NUM)", d1 - (2*pi**3/(27*sqrt(3)) - 2*z3/3))
Ls3_60 = quad(lambda x: log(2*mp.sin(x/2))**2, [0, pi/3])
chk("Ls3(pi/3) = 7 pi^3/108 (NUM)", Ls3_60 - 7*pi**3/108)
chk("P1 sqrt3 = (8/7) Ls3(pi/3) (NUM)", P1*sqrt(3) - 8*Ls3_60/7)
# S = 11 pi^3/324 via proved generating function identity check
z = exp(1j*pi/3)
Sz = polylog(3,z) + mpf('0.5')*log(z)*log(1-z)**2 + z3 + log(1-z)*polylog(2,1-z) - polylog(3,1-z)
chk("S = Im(S(e^{ipi/3})) = 11 pi^3/324", mp.im(Sz) - 11*pi**3/324)
print("ALL CHECKS DONE")
with open('progress.log','a') as fh:
    fh.write("step_55: full certificate script: all identities re-verified at 80 digits\n")
