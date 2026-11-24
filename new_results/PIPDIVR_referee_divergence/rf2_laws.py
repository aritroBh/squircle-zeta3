# rf2: tests of Theorems A(ii), B, C on N_k from rf1 (own code).  Usage: rf2_laws.py NFILE KMAX
import sys, re
from flint import arb, acb, ctx
from rf_series import mul, log, inv
import mpmath as mp
ctx.prec = 600; mp.mp.dps = 120
NF = sys.argv[1]; KM = int(sys.argv[2])

Nmid = {}; Nrad = {}
for line in open(NF):
    k, rest = line.split(" ", 1)
    m = re.match(r"\[(\S+) \+/- (\S+)\]", rest.strip())
    if m:
        Nmid[int(k)] = mp.mpf(m.group(1)); Nrad[int(k)] = mp.mpf(m.group(2))
    else:
        Nmid[int(k)] = mp.mpf(rest.strip()); Nrad[int(k)] = mp.mpf(0)
K = min(KM, max(Nmid) + 1)
pi = mp.pi; g = mp.euler; l2 = mp.log(2)

# ---- closed forms N_0..N_5 (D12) ----
z3, z5 = mp.zeta(3), mp.zeta(5)
cf = [-l2, l2**2/2, pi**2/24 - l2**3/6,
      l2**4/24 + pi**2/24*(1-l2) + z3/4,
      -l2**5/120 + pi**2/24*(1-l2+l2**2/2) + z3/4*(1-l2) + pi**4/2880,
      l2**6/720 + pi**2/24*(1-l2+l2**2/2-l2**3/6) + z3/4*(1-l2+l2**2/2) + pi**4/2880*(1-l2) - 7*pi**4/1920 - z3*pi**2/48 + z5/8]
print("== closed forms N_0..N_5: |N_k(rf1) - closed form| ==")
for k in range(6):
    print(k, mp.nstr(abs(Nmid[k]-cf[k]), 3), " ball radius", mp.nstr(Nrad[k], 3))

# ---- r5 table of the earlier referee (13 printed digits) ----
r5 = "../REF_PIP_referee/r5_log.txt"
print("== vs REF_PIP r5_log (13 digits) ==")
mx = 0
for line in open(r5):
    m = re.match(r"k=\s*(\d+)\s+N_k=\s*(\S+)", line)
    if m:
        k = int(m.group(1)); v = mp.mpf(m.group(2))
        rel = abs(Nmid[k]-v)/abs(v); mx = max(mx, rel)
print("max rel diff k=2..30:", mp.nstr(mx, 3))

# ---- T_k exactly (n = 1 symbol A = e^{2y}(1+e^y)^{-1-h}) ----
Ns = K + 2
fac = [arb(1)]
for n in range(1, 2*Ns+2): fac.append(fac[-1]*n)
ey = [1/fac[n] for n in range(Ns)]
u = [ey[0]+1] + ey[1:]
ell = log(u, Ns, arb(2).log())                  # log(1+e^y)
e2y = [arb(2)**n/fac[n] for n in range(Ns)]
base = mul(e2y, inv(u, Ns), Ns)                  # e^{2y}/(1+e^y)
a = [[None]*Ns for _ in range(Ns)]               # a[m][i] = [y^m h^i] A
P = [arb(1)] + [arb(0)]*(Ns-1)
mell = [-c for c in ell]
for i in range(Ns):
    t = mul(base, P, Ns)
    for m in range(Ns): a[m][i] = t[m]/fac[i]
    P = mul(P, mell, Ns)
def tomp(x):
    return mp.mpf(x.mid().str(110, radius=False))
T = {}
for k in range(2, K):
    s = arb(0)
    for m in range(0, k-1):
        s += fac[m]*a[m][k-2-m]
    T[k] = tomp(s)

# ---- singular approximants T_k^{(L)} (own construction) ----
y0 = mp.mpc(0, pi); Ly0 = mp.log(y0)
# phi_l(h): Phi(u,h) = e^{2u} exp(-(1+h) log((e^u-1)/u)) = sum_l phi_l(h) u^l ; polynomials in h
Lmax = 4; NU = Lmax+1
mp_fac = [mp.factorial(n) for n in range(400)]
phis = [mp.mpf(1)/mp_fac[n+1] for n in range(NU)]          # (e^u-1)/u
# log series via mpmath
def mlog(a, N):
    da = [a[i]*i for i in range(1, len(a))]
    ia = [1/a[0]]
    for n in range(1, N):
        ia.append(-sum(a[k]*ia[n-k] for k in range(1, min(n, len(a)-1)+1))/a[0])
    d = [sum(da[i]*ia[n-i] for i in range(0, n+1) if i < len(da)) for n in range(N-1)]
    return [mp.log(a[0])] + [d[i]/(i+1) for i in range(N-1)]
def mmul(a, b, N):
    return [sum(a[i]*b[n-i] for i in range(0, n+1) if i < len(a) and n-i < len(b)) for n in range(N)]
def mexp(a, N):
    out = [mp.exp(a[0])]
    for n in range(1, N):
        out.append(sum(k*a[k]*out[n-k] for k in range(1, n+1) if k < len(a))/n)
    return out
lphi = mlog(phis, NU)
# Phi = e^{2u} * exp(-lphi) * exp(-h lphi): phi_l(h) = sum_i h^i [u^l] e^{2u}e^{-lphi}(-lphi)^i/i!
c0 = mmul([mp.mpf(2)**n/mp_fac[n] for n in range(NU)], mexp([-c for c in lphi], NU), NU)
phi = [[mp.mpf(0)]*(NU) for _ in range(NU)]   # phi[l][i]
Pw = [mp.mpf(1)] + [mp.mpf(0)]*(NU-1)
for i in range(NU):
    t = mmul(c0, Pw, NU)
    for l in range(NU): phi[l][i] = t[l]/mp_fac[i]
    Pw = mmul(Pw, [-c for c in lphi], NU)
print("phi_0 =", [mp.nstr(x, 6) for x in phi[0][:3]], " phi_1 =", [mp.nstr(x, 6) for x in phi[1][:3]], "(claim phi_1=(3-h)/2)")
def poch_poly(c, m, deg):     # (c+h)_m as coefficient list in h, truncated
    p = [mp.mpf(1)]
    for j in range(m):
        q = [mp.mpf(0)]*(min(len(p)+1, deg+1))
        for i, pi_ in enumerate(p):
            if i < len(q): q[i] += pi_*(c+j)
            if i+1 < len(q): q[i+1] += pi_
        p = q
    return p + [mp.mpf(0)]*(deg+1-len(p))
def TL(k, L):
    tot = mp.mpc(0)
    for l in range(L):
        for m in range(0, k-1):
            i = k-2-m
            ex = [(-Ly0)**n/mp_fac[n] for n in range(i+1)]
            pp = poch_poly(1-l, m, i)
            ph = phi[l][:i+1] + [0]*(i+1-len(phi[l][:i+1]))
            c = mmul(mmul(ph, ex, i+1), pp, i+1)[i]
            tot += (-1)**l * c * y0**(l-1-m)
    return 2*tot.real

def psi(n, x): return mp.psi(n, x)
print("\n== Table 1: Theorem A(ii) split, U_k = N_k - T_k - (-log2)^{k+1}/(k+1)! ==")
print("k  T_k/N_k-1   U_k(2pi)^{k-1}/(k-2)!   U_k/[k(k-2)!6^-k + 2^k]")
for k in list(range(10, K, 10)) + [K-1]:
    U = Nmid[k] - T[k] - (-l2)**(k+1)/mp_fac[k+1]
    print(k, mp.nstr(T[k]/Nmid[k]-1, 3), mp.nstr(U*(2*pi)**(k-1)/mp_fac[k-2], 5), mp.nstr(U/(k*mp_fac[k-2]*mp.mpf(6)**(-k)+mp.mpf(2)**k), 3))

print("\n== Table 2: Theorem C laws ==")
print("k  nu_k  2-term  3-term  r2=(nu-nu2)k^2/log^2k  r3=(nu-nu3)k^3/log^3k   [odd: predicted r2 limit -sigma*pi^2/2]")
rows = []
for k in list(range(11, K)):
    nu = Nmid[k]*pi**(k-1)/(2*mp_fac[k-2])
    # two-term law from the complex formula
    two = 2*(mp_fac[k-2]*y0**(1-k) + mp_fac[k-3]*y0**(2-k)*(psi(0, k-2)+g-Ly0)).real
    nu2 = two*pi**(k-1)/(2*mp_fac[k-2])
    X = lambda j: psi(0, j)+g-mp.log(pi)
    if k % 2:
        sg = (-1)**((k-1)//2)
        nu2b = sg*(1+pi**2/(2*(k-2)))
        nu3 = sg*(1+pi**2/(2*(k-2)) - pi**2*(X(k-3)**2-5*pi**2/12-3+psi(1, k-3))/(2*(k-2)*(k-3)))
    else:
        tau = (-1)**((k-2)//2)
        nu2b = tau*pi*X(k-2)/(k-2)
        nu3 = tau*(pi/(k-2))*(X(k-2)+pi**2*X(k-3)/(2*(k-3)))
    assert abs(nu2-nu2b) < mp.mpf(10)**(-50)
    lk = mp.log(k)
    rows.append((k, nu, nu2, nu3, (nu-nu2)*k**2/lk**2, (nu-nu3)*k**3/lk**3))
for r in rows:
    if r[0] in (11,12,21,22,25,26,29,30,39,40,59,60,61,62,79,80,99,100,119,120,139,140,149,150) or r[0] == K-1:
        print(r[0], mp.nstr(r[1], 10), mp.nstr(r[2], 10), mp.nstr(r[3], 10), mp.nstr(r[4], 4), mp.nstr(r[5], 4))
print("max |r2| over k>=21 odd:", mp.nstr(max(abs(r[4]) for r in rows if r[0] >= 21 and r[0] % 2), 4),
      " even:", mp.nstr(max(abs(r[4]) for r in rows if r[0] >= 21 and r[0] % 2 == 0), 4))
print("max |r3| over k>=21 odd:", mp.nstr(max(abs(r[5]) for r in rows if r[0] >= 21 and r[0] % 2), 4),
      " even:", mp.nstr(max(abs(r[5]) for r in rows if r[0] >= 21 and r[0] % 2 == 0), 4))

print("\n== Table 3: even-k bracket  B_k = tau N_k pi^{k-2}/(2(k-3)!) vs X_{k-2} = psi(k-2)+gamma-log pi ==")
print("k  B_k  X_{k-2}  e=B-X  e*k/log^2k (claim bounded)  e*k/log k  e*(k-3)/X_{k-3} (pred -> pi^2/2=4.93)")
for k in range(12, K, 2):
    tau = (-1)**((k-2)//2)
    Bk = tau*Nmid[k]*pi**(k-2)/(2*mp_fac[k-3]); X2 = psi(0, k-2)+g-mp.log(pi); X3 = psi(0, k-3)+g-mp.log(pi)
    e = Bk - X2
    if k % 10 == 0 or k in (12, 22) or k == K-1 or k == K-2:
        print(k, mp.nstr(Bk, 10), mp.nstr(X2, 10), mp.nstr(e, 6), mp.nstr(e*k/mp.log(k)**2, 4), mp.nstr(e*k/mp.log(k), 4), mp.nstr(e*(k-3)/X3, 5))
print("\n== Table 3b: nu_k / [tau pi (log k + gamma - log pi)/k] - 1, scaled by k/log k (claim bounded) and by k ==")
for k in range(20, K, 20):
    tau = (-1)**((k-2)//2)
    nu = Nmid[k]*pi**(k-1)/(2*mp_fac[k-2]); ref = tau*pi*(mp.log(k)+g-mp.log(pi))/k
    d = nu/ref-1
    print(k, mp.nstr(d, 5), mp.nstr(d*k/mp.log(k), 4), mp.nstr(d*k, 4))

print("\n== Table 4: Theorem B, |T_k - T_k^(L)| / (2(k-2)!/pi^{k-1}) scaled by k^{L+1}/log k (claim bounded) ==")
for k in [21, 22, 41, 42, 61, 62, 81, 82, 101, 102, 121, 122, 141, 142]:
    if k >= K: continue
    sc = 2*mp_fac[k-2]/pi**(k-1)
    out = [k]
    for L in (1, 2, 3, 4):
        d = abs(T[k]-TL(k, L))/sc
        out.append("%s(%s)" % (mp.nstr(d, 3), mp.nstr(d*k**(L+1)/mp.log(k), 3)))
    print(*out)
