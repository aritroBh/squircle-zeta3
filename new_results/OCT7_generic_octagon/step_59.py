code = r'''
import mpmath as mp, sympy as sp, pickle, time, itertools
from fractions import Fraction
mp.mp.dps = 470
t_start = time.time()

# ---------- 1. constants ----------
def Li21(a, c, digits=460):
    a, c = mp.mpf(a), mp.mpf(c)
    N = 2
    while abs(a)**(N+1)*abs(c)**N >= mp.mpf(10)**(-digits-20): N += 1
    apw=[mp.mpf(1)]; cpw=[mp.mpf(1)]
    for i in range(N+2): apw.append(apw[-1]*a); cpw.append(cpw[-1]*c)
    tot = mp.mpf(0)
    for n in range(2, N+2):
        an=apw[n]
        tot += mp.fsum([an*cpw[m]/(mp.mpf(n)*mp.mpf(n)*m) for m in range(1,n)])
    return tot

x1 = Fraction(-1,3); x2 = Fraction(-1,9)
b1 = Fraction(1,4);  b2 = Fraction(1,10)
Lam = {'Lam1': Li21(mp.mpf(x1.numerator)/x1.denominator, mp.mpf((x2/x1).numerator)/(x2/x1).denominator),
       'Lam2': Li21(mp.mpf(x2.numerator)/x2.denominator, mp.mpf((x1/x2).numerator)/(x1/x2).denominator),
       'Lam3': Li21(mp.mpf(b1.numerator)/b1.denominator, mp.mpf((b2/b1).numerator)/(b2/b1).denominator),
       'Lam4': Li21(mp.mpf(b2.numerator)/b2.denominator, mp.mpf((b1/b2).numerator)/(b1/b2).denominator)}
print(f"[{time.time()-t_start:.0f}s] Lam1..4 computed")

# ---------- 2. exact rational point set ----------
def orbit(x):
    return {x, 1-x, 1/x, x/(x-1), 1/(1-x), (x-1)/x}
seeds = [x1, x2, Fraction(-13,27), Fraction(-1,5), Fraction(1,6), Fraction(-7,9), Fraction(-19,81),
         Fraction(-1,2), Fraction(-2,3), Fraction(1,2)]
pts = set()
for s in seeds: pts |= orbit(s)
pts = sorted([p for p in pts if p < 1 and p != 0])
print("Li points (%d):" % len(pts), [str(p) for p in pts])
assert len(pts) == len(set(pts)), "duplicates!"
def f2m(q): return mp.mpf(q.numerator)/q.denominator

l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2 = mp.pi**2/6; z3 = mp.zeta(3)
logs = [(l2,'l2'),(l3,'l3'),(l5,'l5')]
t0=time.time()
Li2v = {p: mp.polylog(2, f2m(p)) for p in pts}
Li3v = {p: mp.polylog(3, f2m(p)) for p in pts}
print(f"[{time.time()-t0:.0f}s] Li2/Li3 at {len(pts)} exact rational points")

basis = [(mp.mpf(1),'1')] + logs + [(z2,'z2')]
for p in pts: basis.append((Li2v[p], f'Li2({p})'))
for (a,na),(b,nb) in itertools.product(logs, repeat=2): basis.append((a*b, na+nb))
basis.append((z3,'z3'))
for p in pts: basis.append((Li3v[p], f'Li3({p})'))
for a,na in logs: basis.append((z2*a, 'z2'+na))
for p in pts:
    for a,na in logs: basis.append((Li2v[p]*a, f'Li2({p}){na}'))
for (a,na),(b,nb),(c,nc) in itertools.product(logs, repeat=3): basis.append((a*b*c, na+nb+nc))
print("basis size:", len(basis))
# near-duplicate guard: all pairwise differences must be > 1e-100
vals = [v for v,_ in basis]
bad = 0
for i in range(len(vals)):
    for j in range(i+1, len(vals)):
        if abs(vals[i]-vals[j]) < mp.mpf(10)**(-100): bad += 1; print("NEAR-DUP:", basis[i][1], basis[j][1])
print("near-duplicate pairs:", bad)
pickle.dump({'basis_names':[n for _,n in basis]}, open('basis_names_exact.pkl','wb'))

# ---------- 3. LLL relation finder with tolerance-stable verification ----------
D = sp.Integer(10)**450
def find_relation(target):
    valsL = vals + [target]
    n = len(valsL)
    rows = []
    for i in range(n):
        r = [sp.Integer(0)]*n + [sp.Integer(int(mp.floor(valsL[i]*D + mp.mpf('0.5'))))]
        r[i] = sp.Integer(1)
        rows.append(r)
    M = sp.Matrix(rows)
    L = M.lll(delta=sp.Rational(3,4))
    Lr = L.tolist()
    rels = []
    for r in Lr:
        if all(x==0 for x in r[:-1]): continue
        if abs(int(r[-1])) < 10**250:
            res = mp.fsum([int(r[i])*valsL[i] for i in range(n)])
            if abs(res) < mp.mpf(10)**(-430):
                rels.append((r, res))
    nonrel_firstn = min((sum(int(x)**2 for x in r[:-1]) for r in Lr if abs(int(r[-1])) >= 10**250), default=None)
    return rels, nonrel_firstn, n

import math
results = {}
for nm, Lv in Lam.items():
    t0=time.time()
    rels, bound, n = find_relation(Lv)
    dt=time.time()-t0
    tgt_rels = [(r,res) for r,res in rels if r[n-1] != 0]
    if tgt_rels:
        r,res = tgt_rels[0]
        # tolerance-stability: verify at 700 digits
        mp.mp.dps = 700
        tgt700 = Li21(*({'Lam1':(x1,x2/x1),'Lam2':(x2,x1/x2),'Lam3':(b1,b2/b1),'Lam4':(b2,b1/b2)}[nm]), )
        # recompute basis at 700
        l2h, l3h, l5h = mp.log(2), mp.log(3), mp.log(5)
        z2h, z3h = mp.pi**2/6, mp.zeta(3)
        lm = {'l2':l2h,'l3':l3h,'l5':l5h}
        import re as _re
        def bval(name):
            if name=='1': return mp.mpf(1)
            if name in lm: return lm[name]
            if name=='z2': return z2h
            if name=='z3': return z3h
            val = mp.mpf(1); rest = name
            if rest.startswith('z2'): val*=z2h; rest=rest[2:]
            mm = _re.match(r'^Li([23])\((-?\d+/\d+|-?\d+)\)', rest)
            if mm:
                q = Fraction(mm.group(2))
                val *= mp.polylog(int(mm.group(1)), f2m(q)); rest = rest[mm.end():]
            for ch in _re.findall(r'l[235]', rest): val *= lm[ch]
            return val
        res700 = mp.fsum([int(r[i])*bval(basis[i][1]) for i in range(n-1)]) + int(r[n-1])*tgt700
        mp.mp.dps = 470
        print(f"{nm}: candidate relation residual@470={mp.nstr(abs(res),3)} residual@700={mp.nstr(abs(res700),3)} [{dt:.0f}s]")
        if abs(res700) < mp.mpf(10)**(-680):
            print(f"  {nm}: REDUCES (tolerance-stable)")
            results[nm] = ('REDUCES', [int(x) for x in r])
            continue
        print(f"  {nm}: candidate SPURIOUS (fails at 700 digits)")
    B = math.isqrt(bound) if bound else None
    print(f"{nm}: SURVIVES at 450dps over basis {n-1}+target (no relation with first-n norm^2 < {bound}, i.e. maxcoeff >~ {B}) [{dt:.0f}s]")
    results[nm] = ('SURVIVES', str(bound))
pickle.dump({'Lam': {k: mp.nstr(v, 460) for k,v in Lam.items()},
             'results': {k:(v[0], str(v[1])) for k,v in results.items()}}, open('lam_lll_450.pkl','wb'))
with open("progress.log","a") as f:
    f.write(f"STEP 06 (exact-rational basis): " + "; ".join(f"{k}:{v[0]}" for k,v in results.items()) + "\n")
print(f"total time {time.time()-t_start:.0f}s")
'''
with open("step_06.py","w") as f: f.write(code)
import subprocess, sys
r = subprocess.run([sys.executable, "step_06.py"], capture_output=True, text=True, cwd="research/08_theory/OCT7_generic_octagon", timeout=3600)
print(r.stdout)
if r.stderr: print("STDERR:", r.stderr[-1000:])
