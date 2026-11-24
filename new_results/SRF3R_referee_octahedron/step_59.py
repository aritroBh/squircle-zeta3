code = r'''"""
verify_srf3r.py — SRF3R referee's consolidated independent verification (< 5 min).
Referee of SRF3_octahedron_radius_proof (claims A-D). All checks use the referee's own code.
 (R1) Theorem A: S(1+d) by DIRECT surface quadrature (graph + radial parametrizations,
      not via T) vs 8(1-tau)^2 T(tau)   [d = +-0.01]
 (R2) zeros of Q: tail min_s |tau*|x/R >= 1 (x>=5); moderate-x min < 1  => S3.2 contour
      claim false as written (V-split repair); corner fiber min |z| >= R, =2R at w=1/2
 (R3) SIGN RESOLUTION (planner finding): raw ratios for the d-series with z0 vs zbar0,
      and C1/C2 fits: Claim D's first d-display sign is wrong; second display correct
 (R4) growth law: beta profile; C_tau, C_d fits; (R,theta) sharpness
 (R5) mixing identity tau->d; exact low-order coefficients vs paper thm:oct;
      3-consecutive-k lemma (S3.4)
 (R6) Pade-Borel resummation vs direct quadrature at d = 0.01
Prints one line per item.
"""
import os, json, time, math
import mpmath as mp
t0 = time.time()
mp.mp.dps = 50
base = os.path.dirname(os.path.abspath(__file__))
srf3 = os.path.join(base, "..", "SRF3_octahedron_radius_proof")
srf1 = os.path.join(base, "..", "SRF1_radii_3d")
z0 = (-mp.log(2) + mp.pi*1j)/2
R = abs(z0); th = mp.arg(z0)
ok = {}

# ---------------- R1: Theorem A by direct quadrature ----------------
def S_graph(p):
    p = mp.mpf(p)
    def integrand(x, y):
        z = (1 - x**p - y**p)**(1/p)
        zx = -x**(p-1)*z**(1-p); zy = -y**(p-1)*z**(1-p)
        return mp.sqrt(1 + zx**2 + zy**2)
    def inner(x):
        ymax = (1 - x**p)**(1/p)
        return ymax * mp.quad(lambda v: integrand(x, ymax*v), [0, 1])
    return 8 * mp.quad(inner, [0, 1])
def S_radial(p):
    p = mp.mpf(p)
    def Nrho(th_, ph):
        st, ct = mp.sin(th_), mp.cos(th_)
        a = (st*mp.cos(ph))**p + (st*mp.sin(ph))**p + ct**p
        rho = a**(-1/p)
        dNth = p*(st*mp.cos(ph))**(p-1)*ct*mp.cos(ph) + p*(st*mp.sin(ph))**(p-1)*ct*mp.sin(ph) - p*ct**(p-1)*st
        dNph = p*(st*mp.cos(ph))**(p-1)*(-st*mp.sin(ph)) + p*(st*mp.sin(ph))**(p-1)*(st*mp.cos(ph))
        drhoth = -(1/p)*a**(-1/p-1)*dNth; drhoph = -(1/p)*a**(-1/p-1)*dNph
        st_ = st if st != 0 else mp.mpf(10)**(-mp.mp.dps)
        return rho*mp.sqrt(rho**2 + drhoth**2 + (drhoph/st_)**2)*st
    return 16 * mp.quad(lambda ph: mp.quad(lambda th_: Nrho(th_, ph), [0, mp.pi/2]), [0, mp.pi/4])
def T_simplex(tau):
    tau = mp.mpf(tau)
    f = lambda t1, t2: ((t1*t2*(1-t1-t2))**(-tau)
                        * mp.sqrt(t1**(2*tau) + t2**(2*tau) + (1-t1-t2)**(2*tau)))
    def inner(t1):
        return (1-t1)*mp.quad(lambda v: f(t1, (1-t1)*v), [0, 1])
    return mp.quad(inner, [0, 1])
worst = 99.0
for ds in ["0.01", "-0.01"]:
    p = 1 + mp.mpf(ds); tau = 1 - 1/p
    St = 8*(1-tau)**2*T_simplex(tau)
    Sg = S_graph(p); Sr = S_radial(p)
    worst = min(worst, float(-mp.log10(abs(Sg-St)/abs(St))), float(-mp.log10(abs(Sr-St)/abs(St))))
ok["R1 Theorem A: graph & radial vs 8(1-tau)^2 T, worst digits (d=+-0.01)"] = round(worst, 1)

# ---------------- R2: zeros of Q / corner fiber ----------------
mp.mp.dps = 35
def Q(t, s, u):
    return (1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + u**(2*t)
def dQ(t, s, u):
    return (2*(1-u)**(2*t)*(s**(2*t)*mp.log(s)+(1-s)**(2*t)*mp.log(1-s))
            + 2*mp.log(1-u)*(1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + 2*mp.log(u)*u**(2*t))
def zero(s, u, m=-1):
    x = -mp.log(u)
    t = (-mp.log(2) + 1j*mp.pi*(2*m+1))/(2*x)
    for _ in range(200):
        tn = t - Q(t, s, u)/dQ(t, s, u)
        if abs(tn-t) < mp.mpf(10)**(-28): return tn
        t = tn
    return t
tail, mid = [], []
for xe, tgt in [(6, tail), (12, tail), (2, mid)]:
    u = mp.e**(-xe); mn = mp.mpf(9)
    for i in range(5, 96, 5):
        s = mp.mpf(i)/100
        for m in [-1, 0]:
            v = abs(zero(s, u, m))*xe/R
            if v < mn: mn = v
    tgt.append(float(mn))
def fib_min(w):
    best = None
    for re0 in [-2,-1,-0.5,0.5,1,2]:
        for im0 in [0.5,1,2,3,5,7]:
            z = mp.mpc(re0, im0)
            for _ in range(100):
                try:
                    F = mp.e**(-2*z*w)+mp.e**(-2*z*(1-w))+1
                    dF = -2*w*mp.e**(-2*z*w)-2*(1-w)*mp.e**(-2*z*(1-w))
                    zn = z - F/dF
                except (OverflowError, ZeroDivisionError):
                    break
                if abs(zn) > 4*R: break
                if abs(zn-z) < mp.mpf(10)**(-25): z = zn; break
                z = zn
            try:
                if abs(mp.e**(-2*z*w)+mp.e**(-2*z*(1-w))+1) < mp.mpf(10)**(-12):
                    if best is None or abs(z) < abs(best): best = z
            except OverflowError:
                pass
    return float(abs(best)/R)
fm = [fib_min(mp.mpf(w)) for w in [0.5, 0.3, 0.1, 0.05, 0.02]]
ok["R2a tail min|tau*|x/R (x=6,12; need>=1)"] = [round(v,4) for v in tail]
ok["R2b moderate-x min|tau*|x/R (x=2; S3.2 contour fails => V-split repair)"] = [round(v,4) for v in mid]
ok["R2c corner fiber min|z|/R (w=0.5,0.3,0.1,0.05,0.02; >=1, =2 at 1/2)"] = [round(v,4) for v in fm]

# ---------------- R3: sign resolution + fits ----------------
mp.mp.dps = 70
j120 = json.load(open(os.path.join(srf3, "SRF3_hi_coeffs_K120.json")))
j60 = json.load(open(os.path.join(srf1, "SRF1_n3_K60.json")))
A_tau = [mp.mpf(x) for x in j120['eta']]
A_d = [mp.mpf(x) for x in j120['d']]
same = all(mp.mpf(j60['d'][k]) == A_d[k] and mp.mpf(j60['tau'][k]) == A_tau[k] for k in range(61))
bt = [A_tau[k]/mp.factorial(k) for k in range(121)]
bd = [A_d[k]/mp.factorial(k) for k in range(121)]
rat = [float((bd[k]/mp.factorial(1))/(mp.mpf(k)**mp.mpf('-1.5')*mp.im(z0**(mp.mpf('0.5')-k)))) for k in range(56, 61)]
ok["R3a SRF1 d,tau == SRF3 d,eta (k<=60)"] = same
ok["R3b d-ratio with Im[z0^(1/2-k)], k=56..60 (theory -sqrt(3pi/4)=%.5f)" % float(-mp.sqrt(3*mp.pi/4))] = [round(v,4) for v in rat]
def linfit(b, beta, k0, k1=120, ncorr=4):
    rows = []; y = []
    for k in range(k0, k1+1):
        c = mp.cos((k-mp.mpf('0.5'))*th); s = mp.sin((k-mp.mpf('0.5'))*th)
        row = [c, s]
        for jj in range(1, ncorr+1): row += [c/k**jj, s/k**jj]
        rows.append(row); y.append(b[k]*R**k*k**(-beta))
    n = len(rows[0]); A = mp.matrix(rows); Y = mp.matrix(y)
    ATA = mp.matrix(n, n); ATY = mp.matrix(n, 1)
    for i in range(n):
        for jj in range(n): ATA[i, jj] = sum(A[r, i]*A[r, jj] for r in range(len(rows)))
        ATY[i] = sum(A[r, i]*Y[r] for r in range(len(rows)))
    x = mp.lu_solve(ATA, ATY)
    res = mp.sqrt(sum((sum(A[r, i]*x[i] for i in range(n))-Y[r])**2 for r in range(len(rows)))/len(rows))
    return [x[i] for i in range(n)], res
xt, rt_ = linfit(bt, mp.mpf('-1.5'), 80)
xd, rd_ = linfit(bd, mp.mpf('-1.5'), 80)
ok["R3c tau-law C1 (theory -sqrt(3 pi R/8) = %s)" % mp.nstr(-mp.sqrt(3*mp.pi*R/8), 14)] = mp.nstr(xt[0], 14)
ok["R3d d-law C2 (theory +sqrt(3 pi R/4) = %s); C1~0" % mp.nstr(mp.sqrt(3*mp.pi*R/4), 14)] = (mp.nstr(xd[1], 14), mp.nstr(xd[0], 4))

# ---------------- R4: beta profile + R/theta sharpness ----------------
prof = []
for ib in [-1.51, -1.505, -1.5, -1.495, -1.49]:
    _, rr = linfit(bt, mp.mpf(ib), 70, 120, 3)
    prof.append((ib, float(rr)))
ok["R4a beta profile (resid: sharp min at -1.5)"] = [(b, "%.1e" % r) for b, r in prof]
def resid_R(iR):
    Rf = R*(1+mp.mpf(iR)/1000000)
    rows = []; y = []
    for k in range(70, 121):
        c = mp.cos((k-mp.mpf('0.5'))*th); s = mp.sin((k-mp.mpf('0.5'))*th)
        row = [c, s, c/k, s/k, c/k**2, s/k**2]
        rows.append(row); y.append(bt[k]*Rf**k*k**mp.mpf('1.5'))
    n = 6; A = mp.matrix(rows); Y = mp.matrix(y)
    ATA = mp.matrix(n, n); ATY = mp.matrix(n, 1)
    for i in range(n):
        for jj in range(n): ATA[i, jj] = sum(A[r, i]*A[r, jj] for r in range(len(rows)))
        ATY[i] = sum(A[r, i]*Y[r] for r in range(len(rows)))
    x = mp.lu_solve(ATA, ATY)
    return float(mp.sqrt(sum((sum(A[r, i]*x[i] for i in range(n))-Y[r])**2 for r in range(len(rows)))/len(rows)))
ok["R4b R profile resid at R*(1+[-1e-5,0,+1e-5]) (A=1/R=%s)" % mp.nstr(1/R, 12)] = ["%.1e" % resid_R(v) for v in [-10, 0, 10]]

# ---------------- R5: mixing, exact low order, 3-consecutive-k ----------------
from math import comb
worst_m = mp.mpf(0)
for k in range(1, 121):
    ssum = mp.mpf(0)
    for jj in range(1, k+1):
        ssum += A_tau[jj]*((-1)**(k-jj))*comb(k-1, jj-1)
    worst_m = max(worst_m, abs(ssum-A_d[k])/abs(A_d[k]))
ok["R5a tau->d mixing identity worst rel err"] = mp.nstr(worst_m, 3)
ex = {2: 1-mp.pi**2/9, 3: 1+mp.pi**2/9-16*mp.zeta(3)/9, 4: 1-2*mp.pi**2/9+32*mp.zeta(3)/9-16*mp.pi**4/405}
ok["R5b d^2,d^3,d^4 vs paper thm:oct (abs diff)"] = [mp.nstr(abs(A_d[k]-ex[k]), 3) for k in (2, 3, 4)]
mp.mp.dps = 30
alpha = mp.pi - th; bound = float(mp.sin(alpha/2))
bad = 0
for k in range(1, 20001):
    if max(abs(mp.cos((kk-mp.mpf('0.5'))*th)) for kk in (k, k+1, k+2)) < bound: bad += 1
ok["R5c 3-consecutive-k lemma violations (k<=20000)"] = bad

# ---------------- R6: Pade-Borel at d=0.01 ----------------
mp.mp.dps = 50
b = [A_d[k]/mp.factorial(k) for k in range(121)]
def pade_borel_sum(b, d, L, M):
    A = mp.matrix(M, M); Y = mp.matrix(M, 1)
    for i in range(M):
        for jj in range(M):
            idx = L - jj + i
            A[i, jj] = b[idx] if idx >= 0 else mp.mpf(0)
        Y[i] = -b[L + i + 1]
    qv = mp.lu_solve(A, Y)
    q = [mp.mpf(1)] + [qv[i] for i in range(M)]
    p = []
    for k in range(L+1):
        s = mp.mpf(0)
        for jj in range(min(k, M)+1): s += q[jj]*b[k-jj]
        p.append(s)
    def P(z):
        return sum(p[k]*z**k for k in range(L+1))/sum(q[jj]*z**jj for jj in range(M+1))
    return mp.quad(lambda t: mp.e**(-t)*P(t*d), [0, mp.inf])
tau = 1 - 1/mp.mpf('1.01')
St = 8*(1-tau)**2*T_simplex(tau)/(4*mp.sqrt(3))
v = pade_borel_sum(b, mp.mpf('0.01'), 45, 45)
ok["R6 Pade-Borel[45/45] at d=0.01 vs T-based value (digits)"] = round(float(-mp.log10(abs(v-St)/abs(St))), 1)

print("=" * 78)
for k, v in ok.items(): print(f"  {k}: {v}")
print(f"verify_srf3r total time: {time.time()-t0:.1f}s")
'''
open("verify_srf3r.py","w").write(code)
print("written", len(code))
