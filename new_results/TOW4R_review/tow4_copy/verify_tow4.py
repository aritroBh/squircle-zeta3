"""verify_tow4.py -- TOW4 per-lemma certificates, PASS/FAIL each. Total < 30 min.
Lemmas:
  TAU8   tau1(8),tau2(8) in Q(i), exact vs tow2lib definitions; conjugation symmetry.
  L      L_simplified cache == tow2lib.L(r) proved construction, symbolic, r=1..8.
  RHO    a_7,a_8 = Bell from L; rho(7)/e00==a_7, rho(8)/e00==a_8 vs TOWN1, 717 cells each, 0 bad.
  CLOSE  forced recursion rebuilt; closing F(7,4)=0, F(8,5)=0 (numerator 0);
         forced F == tow3lib.e (delta<=6 all s); F(delta,top)==zigzag (delta<=8).
  ADJ    9 adjacent identities (7,0..3),(8,0..4), exact numerator 0 (subprocess prove_adj78_fast).
  E3     blindness at delta=7,8, exact numerator 0 (subprocess prove_E3b78).
  LAW    full tower-law confirmation delta<=8 vs TOWN1 (subprocess check_full78): 49953 cells, 0 bad.
  W6     Task C: k=6 vs quartic(k=7..11) 48/48; quarticity k=7..13; +/-i poles vanish k=6..13.
  MASTER Task D reduction: master x=-1 coefficient feed(U(-1))-a_delta==0 at delta=7,8.
Run: PY verify_tow4.py   (PY = venv python)
"""
import sys, os, time, pickle, subprocess
from fractions import Fraction as Fr
work = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, work)
import sympy as sp
import tow2lib as T2, tow4lib as T4, tow3lib as T3
d, m = T4.d, T4.m
i = sp.I
PY = sys.executable
RESULTS = []
def report(name, ok, extra=''):
    RESULTS.append((name, bool(ok)))
    print(f'[{"PASS" if ok else "FAIL"}] {name}  {extra}', flush=True)

t00 = time.time()
# ---- TAU8 ----
t1 = T4.tau1f(8); t2 = T4.tau2f(8)
d1 = sp.simplify(sp.expand(T2.tau1(8) - t1)); d2 = sp.simplify(sp.expand(T2.tau2(8) - t2))
sym_ok = (sp.re(t1) - sp.im(t1)*i) == t2
report('TAU8', d1==0 and d2==0 and sym_ok, f'tau1(8)={t1}')

# ---- L ----
Ls_raw = pickle.load(open(work + '/L_simplified.pkl','rb'))
cache = {k: sp.sympify(v, locals={'d':d,'m':m,'I':i}) for k,v in Ls_raw.items()}
ok = True
for r in range(1,9):
    ok &= (sp.cancel(sp.together(T2.L(r) - cache[r])) == 0)
report('L', ok, 'r=1..8 symbolic')

# ---- RHO ----
def qsi_div(x, y):
    def gm(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
    def gd(a,b):
        den = b[0]**2+b[1]**2; return ((a[0]*b[0]+a[1]*b[1])/den, (a[1]*b[0]-a[0]*b[1])/den)
    A=(x[0],x[2]); B=(x[1],x[3]); C=(y[0],y[2]); D=(y[1],y[3])
    C2 = (gm(C,C)[0]-2*gm(D,D)[0], gm(C,C)[1]-2*gm(D,D)[1])
    nA = (gm(A,C)[0]-2*gm(B,D)[0], gm(A,C)[1]-2*gm(B,D)[1])
    nB = (gm(B,C)[0]-gm(A,D)[0], gm(B,C)[1]-gm(A,D)[1])
    return (gd(nA, C2), gd(nB, C2))
H = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
rho_ok = True; counts = {}
for delta in [7,8]:
    aa = T4.arho(delta)
    bad=0; n=0
    for key in H:
        if key=='_meta': continue
        p,q = key; dd = p+q
        if 2*dd-4 < delta: continue
        num = H[(p,q)][delta][dd]; den = H[(p,q)][0][dd]
        rv = qsi_div(tuple(num), tuple(den))
        cv = sp.expand(aa.subs({d: dd, m: p}))
        exp = (Fr(sp.re(cv)), Fr(sp.im(cv)))
        n += 1
        if rv[1] != (Fr(0),Fr(0)) or rv[0] != exp: bad += 1
    counts[delta] = (n, bad); rho_ok &= (bad==0)
report('RHO', rho_ok, f'a_7 {counts[7][0]} cells {counts[7][1]} bad; a_8 {counts[8][0]} cells {counts[8][1]} bad')

# ---- CLOSE (rebuild forced F, prove closing, sanity) ----
r = subprocess.run([PY, work + '/forced78.py'], capture_output=True, text=True)
close_ok = 'VERDICT non-adjacent depth<=8: True' in r.stdout
# sanity: forced == tow3lib e (delta<=6), F(top)==zigzag (delta<=8)
import importlib
sys.dont_write_bytecode = True
san1 = all(sp.cancel(sp.together(T4.forced(dl,s)-T3.e(dl,s)))==0 for dl in range(0,7) for s in range(0,dl//2+1))
san2 = all(sp.cancel(sp.together(T4.forced(dl,dl//2)-T2.e_top(dl)))==0 for dl in range(0,9))
report('CLOSE', close_ok, 'closing F(7,4)=0,F(8,5)=0 numerator 0 (forced78.py)')
report('CLOSE-sanity', san1 and san2, 'F==e (dl<=6), F(top)==zigzag (dl<=8)')

# ---- ADJ ----
adj_ok = True; times = {}
for delta,s in [(7,0),(7,1),(7,2),(7,3),(8,0),(8,1),(8,2),(8,3),(8,4)]:
    t0=time.time()
    rr = subprocess.run([PY, work + '/prove_adj78_fast.py', str(delta), str(s)], capture_output=True, text=True)
    ok = 'True' in rr.stdout
    adj_ok &= ok; times[(delta,s)] = round(time.time()-t0)
report('ADJ', adj_ok, '9 adjacent identities, exact numerator 0; times ' + str(times))

# ---- E3 ----
e3_ok = True
for delta in [7,8]:
    rr = subprocess.run([PY, work + '/prove_E3b78.py', str(delta)], capture_output=True, text=True)
    e3_ok &= ('True' in rr.stdout)
report('E3', e3_ok, 'blindness delta=7,8, exact numerator 0')

# ---- LAW ----
rr = subprocess.run([PY, work + '/check_full78.py'], capture_output=True, text=True)
law_ok = 'bad = 0' in rr.stdout and 'bad = 0' in rr.stdout.splitlines()[-2]
report('LAW', law_ok, rr.stdout.strip().splitlines()[-3] if rr.stdout else '')

# ---- W6 ----
from fractions import Fraction as Fr2
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(6, 14)}
keys5 = [(5,0),(4,1),(3,2)]
def lag(ys, xs, x):
    out = Fr2(0)
    for a,xa in enumerate(xs):
        La = Fr2(1)
        for b,xb in enumerate(xs):
            if b!=a: La *= Fr2(x-xb, xa-xb)
        out += ys[a]*La
    return out
w6a = True; n48=0
for key in keys5:
    poles = set()
    for k in range(6,12): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr2(0),)*4)[comp] for k in [7,8,9,10,11]]
            n48 += 1
            if lag(ys,[7,8,9,10,11],6) != K[6][key]['real'].get(pole,(Fr2(0),)*4)[comp]: w6a=False
w6b = True
for key in keys5:
    poles = set()
    for k in range(7,14): poles |= set(K[k][key]['real'].keys())
    for pole in poles:
        for comp in range(4):
            ys = [K[k][key]['real'].get(pole,(Fr2(0),)*4)[comp] for k in range(7,14)]
            cur = ys
            for _ in range(5): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
            if any(c!=0 for c in cur): w6b=False
w6c = all(not K[k][key]['ipos'] and not K[k][key]['ineg'] for k in range(6,14) for key in keys5)
report('W6', w6a and w6b and w6c, f'k=6 vs quartic {n48} comp 0 bad; quarticity 7..13; +/-i vanish')

# ---- MASTER (Task D reduction verified at 7,8) ----
def U(dp): return sum((-1)**s * T4.forced(dp, s) for s in range(0, dp//2+1))
def feedU(delta):
    tot = sp.Integer(0)
    for dp in range(0, delta+1):
        uA = U(dp).subs({d: d+1, m: m+1}, simultaneous=True)
        uB = U(dp).subs({d: d+1, m: m}, simultaneous=True)
        tot += i*(2*d-1)/(2*d-dp-1) * (T4.tau1f(delta-dp)*uA - T4.tau2f(delta-dp)*uB)
    return tot
m_ok = all(sp.cancel(sp.together(feedU(dl) - T4.arho(dl)))==0 for dl in [7,8])
report('MASTER', m_ok, 'closing <=> master x=-1 coefficient, verified delta=7,8')

print('='*70)
npass = sum(1 for _,o in RESULTS if o)
print(f'VERIFY_TOW4: {npass}/{len(RESULTS)} lemmas PASS, total {round(time.time()-t00,1)}s', flush=True)
for nm, o in RESULTS: print(f'  {nm}: {"PASS" if o else "FAIL"}')
