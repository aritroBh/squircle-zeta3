verify = r'''"""verify_tow5.py -- TOW5 per-lemma certificates, PASS/FAIL each. Total < 40 min.
Lemmas:
  TAU910   tau1/tau2(9),(10) in Q(i), exact vs tow2lib definitions; conjugation symmetry.
  R5       R_5 = N_5/(3 Pi_5) proved (subprocess prove_S5.py: P_5 + symbolic identity + direct check).
  L        L cache == tow2lib-style L(r), symbolic, r = 1..10 (subprocess verify_L.py).
  RHO      a_9,a_10 = Bell from L; rho(9)/e00==a_9, rho(10)/e00==a_10 vs step-B data, 425 cells each.
  CLOSE    closing F(9,5)=0, F(10,6)=0 rebuilt inline from arho + PROVED delta<=8 family (numerator 0);
           sanity F==tow3lib.e (delta<=6), F(top)==zigzag (delta<=10).
  ADJ      11 adjacent identities (9,0..4),(10,0..5), exact numerator 0 (subprocess prove_adj910_fast).
  E3       blindness at delta=9,10, exact numerator 0 (subprocess prove_E3b910).
  LAW      full tower-law confirmation delta<=10 vs step-B data (subprocess check_full910): 990+11 cells.
  W7       Task D: k=7 vs quintic(8..13) 56/56; quinticity 7..13; +-i vanish; k=6 off (subprocess).
  MASTER   Task E: master x=-1 coefficient feed(U(-1))-a_delta==0 at delta=9,10.
Run: PY verify_tow5.py   (PY = venv python). Subprocess return codes checked (TOW4R erratum).
"""
import sys, os, time, pickle, subprocess
from fractions import Fraction as Fr
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work)
import sympy as sp
import tow2lib as T2, tow5lib as T5, tow3lib as T3
d, m = T5.d, T5.m
i = sp.I
PY = sys.executable
RESULTS = []
def report(name, ok, extra=''):
    RESULTS.append((name, bool(ok)))
    print(f'[{"PASS" if ok else "FAIL"}] {name}  {extra}', flush=True)
def run(name, script, *args, timeout=1700):
    t0 = time.time()
    r = subprocess.run(['nice', '-n', '10', PY, os.path.join(work, script), *map(str, args)],
                       capture_output=True, text=True, timeout=timeout)
    ok = (r.returncode == 0)
    print(f'    [{script} {" ".join(map(str,args))}] rc={r.returncode} ({round(time.time()-t0,1)}s)', flush=True)
    if not ok:
        print('    STDOUT tail:', r.stdout[-400:]); print('    STDERR tail:', r.stderr[-400:])
    return ok

t00 = time.time()
# ---- TAU910 ----
tau910 = pickle.load(open(os.path.join(work, 'tau910.pkl'), 'rb'))
ok = True; msgs = []
for n in [9, 10]:
    t1 = sp.sympify(tau910[n][0], locals={'I': i}); t2 = sp.sympify(tau910[n][1], locals={'I': i})
    w = sp.Symbol('w')
    def qi_diff(e):
        g = sp.together(e); num, den = sp.fraction(g)
        def par(x):
            P = sp.Poly(sp.expand(x.subs({sp.sqrt(2): w})), w); a = sp.Integer(0); b = sp.Integer(0)
            for (e2,), c in P.terms():
                if e2 % 2 == 0: a += c * 2**(e2//2)
                else: b += c * 2**((e2-1)//2)
            return a, b
        na, nb = par(num); da, db = par(den)
        return sp.expand(nb*da - na*db), sp.expand(na*da - 2*nb*db), sp.expand(da**2 - 2*db**2)
    s1 = qi_diff(T2.tau1(n) - t1); s2 = qi_diff(T2.tau2(n) - t2)
    ok &= (s1[0] == 0 and s1[1] == 0 and s2[0] == 0 and s2[1] == 0)
    ok &= (sp.expand(sp.conjugate(t1) - (-1)**n * t2) == 0)
    msgs.append(f'tau1({n})={t1}')
report('TAU910', ok, '; '.join(msgs))

# ---- R5 (subprocess) ----
report('R5', run('R5', 'prove_S5.py'), 'R_5 = N_5/(3 Pi_5); P_5 proved; direct-diff check (prove_S5.py)')

# ---- L (subprocess) ----
report('L', run('L', 'verify_L.py'), 'L cache == tow2lib-style L(r), r = 1..10, symbolic')

# ---- RHO ----
def qsi_div(x, y):
    def gm(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    def gd(a, b):
        den = b[0]**2+b[1]**2; return ((a[0]*b[0]+a[1]*b[1])/den, (a[1]*b[0]-a[0]*b[1])/den)
    A = (x[0], x[2]); B = (x[1], x[3]); C = (y[0], y[2]); D = (y[1], y[3])
    C2 = (gm(C, C)[0]-2*gm(D, D)[0], gm(C, C)[1]-2*gm(D, D)[1])
    nA = (gm(A, C)[0]-2*gm(B, D)[0], gm(A, C)[1]-2*gm(B, D)[1])
    nB = (gm(B, C)[0]-gm(A, D)[0], gm(B, C)[1]-gm(A, D)[1])
    return (gd(nA, C2), gd(nB, C2))
RHO = pickle.load(open(os.path.join(work, 'data_rho_D10_d7_40.pkl'), 'rb'))
rho_ok = True; counts = {}
for delta in [9, 10]:
    aa = T5.arho(delta)
    bad = 0; n = 0
    for dd in RHO:
        if 2*dd-4 < delta: continue
        for key in RHO[dd]:
            p, q = key
            rv = qsi_div(RHO[dd][key][delta], RHO[dd][key][0])
            cv = sp.expand(aa.subs({d: dd, m: p}))
            exp = (Fr(sp.re(cv)), Fr(sp.im(cv)))
            n += 1
            if rv[1] != (Fr(0), Fr(0)) or rv[0] != exp: bad += 1
    counts[delta] = (n, bad); rho_ok &= (bad == 0)
report('RHO', rho_ok, f'a_9 {counts[9][0]} cells {counts[9][1]} bad; a_10 {counts[10][0]} cells {counts[10][1]} bad')

# ---- CLOSE (rebuild closing inline from arho + PROVED delta<=8 family) ----
raw78 = pickle.load(open(os.path.join(work, 'forced_F78.pkl'), 'rb'))
Fz = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': i}) for k, v in raw78.items()}
def F(dp, s):
    if s > dp//2: return sp.Integer(0)
    return Fz[(dp, s)]
def Lam(delta, s):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        fA = F(dp, s).subs({d: d+1, m: m+1}, simultaneous=True)
        fB = F(dp, s).subs({d: d+1, m: m}, simultaneous=True)
        tot += i*(2*d-1)/(2*d-dp-1) * (T5.tau1f(delta-dp)*fA - T5.tau2f(delta-dp)*fB)
    return tot
close_ok = True
for delta in [9, 10]:
    Fz[(delta, 0)] = sp.cancel(T5.arho(delta))
    for s in range(0, delta//2+1):
        nxt = sp.together(Lam(delta, s) - F(delta, s))
        if s+1 <= delta//2:
            Fz[(delta, s+1)] = sp.cancel(nxt)
        else:
            close_ok &= (sp.expand(sp.numer(nxt)) == 0)
san1 = all(sp.cancel(sp.together(T5.forced(dl, s)-T3.e(dl, s))) == 0 for dl in range(0, 7) for s in range(0, dl//2+1))
san2 = all(sp.cancel(sp.together(T5.forced(dl, dl//2)-T2.e_top(dl))) == 0 for dl in range(0, 11))
report('CLOSE', close_ok, 'closing F(9,5)=0, F(10,6)=0 rebuilt inline, numerator 0')
report('CLOSE-sanity', san1 and san2, 'F==e (dl<=6), F(top)==zigzag (dl<=10)')

# ---- ADJ (subprocess, 11 identities) ----
adj_ok = True
for delta, s in [(9,0),(9,1),(9,2),(9,3),(9,4),(10,0),(10,1),(10,2),(10,3),(10,4),(10,5)]:
    adj_ok &= run(f'ADJ({delta},{s})', 'prove_adj910_fast.py', delta, s)
report('ADJ', adj_ok, '11 adjacent identities, exact numerator 0')

# ---- E3 (subprocess) ----
e3_ok = run('E3(9)', 'prove_E3b910.py', 9) & run('E3(10)', 'prove_E3b910.py', 10)
report('E3', e3_ok, 'blindness delta=9,10, exact numerator 0')

# ---- LAW (subprocess) ----
report('LAW', run('LAW', 'check_full910.py'), 'full tower law delta<=10 vs step-B data: 990+11 cells')

# ---- W7 (subprocess) ----
report('W7', run('W7', 'verify_W7.py'), 'k=7==quintic 56/56; quinticity 7..13; +-i vanish; k=6 off 18/56')

# ---- MASTER (Task E verified at 9,10) ----
def U(dp): return sum((-1)**s * T5.forced(dp, s) for s in range(0, dp//2+1))
def feedU(delta):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        uA = U(dp).subs({d: d+1, m: m+1}, simultaneous=True)
        uB = U(dp).subs({d: d+1, m: m}, simultaneous=True)
        tot += i*(2*d-1)/(2*d-dp-1) * (T5.tau1f(delta-dp)*uA - T5.tau2f(delta-dp)*uB)
    return tot
m_ok = all(sp.cancel(sp.together(feedU(dl) - T5.arho(dl))) == 0 for dl in [9, 10])
report('MASTER', m_ok, 'closing <=> master x=-1 coefficient, verified delta=9,10')

print('=' * 70)
npass = sum(1 for _, o in RESULTS if o)
print(f'VERIFY_TOW5: {npass}/{len(RESULTS)} lemmas PASS, total {round(time.time()-t00,1)}s', flush=True)
for nm, o in RESULTS: print(f'  {nm}: {"PASS" if o else "FAIL"}')
'''
work = '<project>/research/08_theory/TOW5_tower_depth10'
open(work+'/verify_tow5.py','w').write(verify)
print("verify_tow5.py written:", len(verify), "chars")
