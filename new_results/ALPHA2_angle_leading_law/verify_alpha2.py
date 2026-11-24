"""verify_alpha2.py -- PASS/FAIL runner for ALPHA2 (corner-angle leading law).
Parts:
  A1: top table = single factor F_d(c) x squircle pattern; F_d t-polynomial fits (8 angles, d<=8)
  S1: symbolic corner-measure identities (sympy, PROVED-level)
  C1: corner Master Lemma (i)(ii)(iii) by quadrature (3 angles x 3 strip points, 40 digits)
  C2: E_alpha closed form = integral definition; squircle reduction; Jfull Legendre form
  B1: the leading law c(w;alpha) at c=1/2 and c=2/3, w=2..8 (40+ digits)
  D1: cross-checks e_2 = 2T/(1-ca), e_3 = -J/(2c^2), e_1 = -ln sin(a/2)/cos^2(a/2) (3 angles)
Expected runtime ~6-9 min. Prints PASS/FAIL per part with worst residuals.
"""
import sys, os, pickle, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from math import factorial
from mpmath import mp, mpf, nstr, pi, sin, cos, acos, gamma, hyper, diffs, quad, sqrt, polylog, zeta, log, legenp

t00 = time.time()
RESULTS = []
def report(name, ok, detail):
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}", flush=True)

def qpart(t): return Fr(*t[0])
def fdiff(seq, n):
    for _ in range(n): seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
    return seq

def extract_beta18(name, dmax=8):
    data = pickle.load(open(f'kern18_{name}.pkl','rb'))
    K = data['kernels']; r = Fr(*data['r']); c = Fr(*data['c'])
    tab = {}
    for d in range(1, dmax+1):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            seqs = {}
            for k in range(d, 19):
                P = K[k][key]
                for rho, cf in P['real'].items(): seqs.setdefault(rho, []).append(qpart(cf))
            for rho, s in seqs.items():
                if len(s) < d + 3: continue
                top = fdiff(list(s[1:]), d-1)
                assert all(x == top[0] for x in top)
                beta = top[0]/factorial(d-1)
                if beta != 0: tab[(d, key, Fr(*rho))] = beta
    return c, r, tab

# ---------- A1: single factor + fits ----------
import sympy as sp
bs = pickle.load(open('beta_squircle.pkl','rb'))
sq = {}
for (d,ij,p),v in bs.items(): sq[(d,ij,{'0':0,'1':1,'-1':-1,'-r':2}[p])] = v[0]
names = ['c12','c35','c13','c45','c23','c57','c14','c67']
nfail = 0; Fd = {}
for name in names:
    c, r, tab = extract_beta18(name)
    for d in range(1, 9):
        ratios = set()
        for (dd,ij,p),v in tab.items():
            if dd!=d: continue
            marker = 2 if p == -r else int(p)
            if (d,ij,marker) in sq and sq[(d,ij,marker)]!=0: ratios.add(v/sq[(d,ij,marker)])
        if len(ratios) > 1: nfail += 1
        elif len(ratios)==1: Fd.setdefault(d,{})[c] = ratios.pop()
report("A1a single factor F_d per (c,d)", nfail==0, f"{nfail} multi-ratio failures, 8 angles d<=8")
# fits
tsym = sp.symbols('t')
expect = {2: tsym, 3: (1+3*tsym)/4, 4: tsym*(2+3*tsym)/5, 5: (2+15*tsym+15*tsym**2)/32,
          6: tsym*(17+60*tsym+45*tsym**2)/122, 7: (17+231*tsym+525*tsym**2+315*tsym**3)/1088,
          8: tsym*(62+378*tsym+630*tsym**2+315*tsym**3)/1385}
bad = 0
for d in range(2, 9):
    pts = [(sp.Rational(1), sp.Rational(1))]
    for c, v in Fd[d].items():
        x = sp.Rational(c.numerator, c.denominator)
        pts.append(((1-x**2)/x**2, sp.Rational(v.numerator, v.denominator)))
    for tp, yp in pts:
        if sp.simplify(expect[d].subs(tsym, tp) - yp) != 0: bad += 1
report("A1b F_d(t) closed forms d<=8", bad==0, f"{bad} point failures over {sum(len(Fd[d])+1 for d in range(2,9))} points")

# ---------- S1: symbolic measures ----------
u = sp.symbols('u'); k = sp.symbols('kappa', positive=True)
S = sp.sqrt(u**2 + 4*k**2); D = sp.sqrt(1+4*k**2)
z = (u+S)/(2*k); gamma_ = (D+1)/(2*k); sroot = (D-1)/(2*k)
mu = {1: sp.Rational(1,2)*(1/S + 1/(u-1) + D/(S*(u-1))),
      -1: sp.Rational(1,2)*(1/S + 1/(u+1) - D/(S*(u+1))),
      -2: sp.Rational(1,2)*(1/S + 1/(u-1) - D/(S*(u-1))), 0: 1/S}
dz_du = sp.diff(z, u)
okm = all(sp.simplify(dz_du/(z - z0) - muu) == 0 for z0, muu in [(gamma_,mu[1]),(-gamma_,mu[-1]),(-sroot,mu[-2]),(0,mu[0])])
m_odd = mu[-1]-mu[1]; m_e = mu[-1]-mu[-2]; m_diag = mu[-1]-mu[0]/2
t1 = sp.simplify(m_odd - (1 + D*u/S)/(1-u**2))
t2 = sp.simplify(m_e - (1 - D/S)/(1-u**2))
t3 = sp.simplify(m_odd - (1/(1-u) - u*m_e))
t4 = sp.simplify(m_diag - (1-u)/2*m_e)
ok2 = all(sp.radsimp(sp.simplify(xx))==0 for xx in [t1,t2,t3,t4])
report("S1 corner measures + LEAD2 identities", okm and ok2, "sympy residuals zero" if okm and ok2 else "nonzero residual")

# ---------- C1: Master Lemma ----------
mp.dps = 45
def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Eint(s, ca): return quad(lambda w: (w**(1+s) + w**(-s))/Qf(w,ca)**3, [0, 1])
def WWq(s, ca): return quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [0, 1]) + quad(lambda w: w**(-s)*(1+w)/Qf(w,ca)**3, [1, mpf('inf')])
def FFq(s, ca):
    h = lambda w: ((1+w)/Qf(w,ca) - 1)/w
    return quad(lambda w: w**(-s)*h(w), [0, 1]) + quad(lambda w: w**(-s)*h(w), [1, mpf('inf')])
def GGq(s, ca):
    g1 = lambda w: (1-w)/(1+w)*((1+w)/Qf(w,ca) - 1)/w
    return quad(lambda w: w**(-s)*g1(w), [0, 1]) + quad(lambda w: w**(-s)*g1(w), [1, mpf('inf')])
worst = {1: mpf(0), 2: mpf(0), 3: mpf(0)}
for alpha in [pi/3, 2*pi/3, pi/6]:
    ca = cos(alpha)
    for s in [mpf('0.3'), mpf('-0.2'), mpf('0.55')]:
        worst[1] = max(worst[1], abs(WWq(s,ca) - (Eint(s,ca)+Eint(-s,ca))))
        worst[2] = max(worst[2], abs(FFq(s,ca) - (1+ca)*(Eint(s,ca)-Eint(-s,ca))/s))
        worst[3] = max(worst[3], abs(GGq(s,ca) - (2*pi/sin(pi*s) - (1-ca)*WWq(s,ca)/s)))
report("C1 Master Lemma (i)(ii)(iii)", all(worst[i] < mpf('1e-18') for i in [1,2,3]),
       f"worst residuals i={nstr(worst[1],2)}, ii={nstr(worst[2],2)}, iii={nstr(worst[3],2)}")

# ---------- C2: E_alpha closed form, squircle reduction, Jfull ----------
def Ecf(s, ca):
    return gamma(1-s)*gamma(2+s)*hyper([1+s, -s],[2], (1+ca)/2)/(1-ca)
w2 = mpf(0)
for alpha in [pi/3, 2*pi/3, pi/6]:
    ca = cos(alpha)
    for s in [mpf('0.3'), mpf('-0.2')]:
        w2 = max(w2, abs(Eint(s,ca) - Ecf(s,ca)))
# squircle reduction
def Esq(s): return gamma(1+s/2)*gamma(mpf('0.5')-s/2)/sqrt(pi)
esq = max(abs(Esq(s) - Ecf(s, mpf('0'))) for s in [mpf('0.3'), mpf('-0.2')])
# Jfull Legendre
def Jtrig(s, alpha): return quad(lambda uu: sin(uu)**(s-1)*sin(alpha+uu)**(-s), [0, pi-alpha])
wj = mpf(0)
for alpha in [pi/3, 2*pi/3]:
    ca = cos(alpha)
    for s in [mpf('0.3'), mpf('0.5')]:
        wj = max(wj, abs(Jtrig(s,alpha) - pi*legenp(s-1, 0, -ca, type=2)/sin(pi*s)))
report("C2 E_alpha closed form + squircle reduction + Jfull", w2 < mpf('1e-30') and esq < mpf('1e-40') and wj < mpf('1e-12'),
       f"Ecf-int={nstr(w2,2)}, squircle={nstr(esq,2)}, Jfull={nstr(wj,2)}")

# ---------- B1: the law at c=1/2, c=2/3, w=2..8 ----------
from assemble import N_integrals
A = [1,1,1,2,5,16,61,272,1385,7936]
def eps(m): return (-1)**((m+1)//2)
def ew_cf(ca, nmax=9):
    return [list(diffs(lambda s: Ecf(s, ca), 0, n))[-1]/factorial(n) for n in range(nmax)]
fits = {2: tsym, 3: (1+3*tsym)/4, 4: tsym*(2+3*tsym)/5, 5: (2+15*tsym+15*tsym**2)/32,
        6: tsym*(17+60*tsym+45*tsym**2)/122, 7: (17+231*tsym+525*tsym**2+315*tsym**3)/1088}
wl = mpf(0); detail = []
for name in ['c12','c23']:
    c, r, tab = extract_beta18(name, dmax=7)
    N = N_integrals(r, dmax=7, dps=45)
    ca = 2*(mpf(c.numerator)/c.denominator)**2 - 1
    ew = ew_cf(ca)
    tv = (1-ca)/(1+ca)
    for w in range(2, 9):
        d = w-1
        tot = mpf(0)
        for (dd, key, rho), beta in tab.items():
            if dd != d: continue
            rhom = mpf(rho.numerator)/rho.denominator
            tot += mpf(beta.numerator)/beta.denominator * N[(key, rhom)]
        if d == 1: Fdv = mpf(1)
        else: Fdv = mpf(str(sp.N(fits[d].subs(tsym, tv), 45)))
        trig = (1-ca)/2 if w % 2 == 0 else (1+ca)/2
        pred = eps(w-2)*A[w-1]*Fdv*trig*ew[w]/factorial(w-2)
        wl = max(wl, abs(tot-pred))
    detail.append(f"{name}: ok")
report("B1 leading law c(w;a), w=2..8, c=1/2,2/3", wl < mpf('1e-38'), f"worst |asm-pred|={nstr(wl,2)}")

# ---------- D1: cross-checks ----------
def Jc(c):
    return log(1-c**2)*(polylog(2,c**2)-zeta(2)) + 2*polylog(3,1-c**2) - 2*zeta(3) + log(c**2)*log(1-c**2)**2
wx = mpf(0)
for c in ['1/2','3/5','2/3']:
    cm = mpf(c); ca = 2*cm**2-1; alpha = acos(ca)
    ew = ew_cf(ca, 4)
    T = pi**2/12 - polylog(2, cm**2)/2
    wx = max(wx, abs(ew[2] - T/((1-ca)/2)))
    wx = max(wx, abs(ew[3] + Jc(cm)/(2*cm**2)))
    wx = max(wx, abs(ew[1] + log(sin(alpha/2))/cos(alpha/2)**2))
report("D1 e_2=2T/(1-ca), e_3=-J/(2c^2), e_1=-ln sin(a/2)/cos^2(a/2)", wx < mpf('1e-40'), f"worst={nstr(wx,2)}")

print(f"\nTOTAL: {sum(1 for _,ok,_ in RESULTS if ok)}/{len(RESULTS)} PASS, {time.time()-t00:.0f}s")
for name, ok, detail in RESULTS:
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
if not all(ok for _, ok, _ in RESULTS):
    sys.exit(1)
