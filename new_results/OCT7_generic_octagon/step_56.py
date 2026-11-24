code = r'''
import mpmath as mp, sympy as sp, pickle, time, itertools, json
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
def Li12(a, c, digits=460):
    a, c = mp.mpf(a), mp.mpf(c)
    N = 2
    while abs(a)**(N+1)*abs(c)**N >= mp.mpf(10)**(-digits-20): N += 1
    apw=[mp.mpf(1)]; cpw=[mp.mpf(1)]
    for i in range(N+2): apw.append(apw[-1]*a); cpw.append(cpw[-1]*c)
    tot = mp.mpf(0)
    for n in range(2, N+2):
        an=apw[n]
        tot += mp.fsum([an*cpw[m]/(n*mp.mpf(m)*mp.mpf(m)) for m in range(1,n)])
    return tot

x1 = mp.mpf(-1)/3; x2 = mp.mpf(-1)/9
b1 = mp.mpf(1)/4;  b2 = mp.mpf(1)/10
t0=time.time()
Lam = {'Lam1': Li21(x1, x2/x1), 'Lam2': Li21(x2, x1/x2), 'Lam3': Li21(b1, b2/b1), 'Lam4': Li21(b2, b1/b2)}
print(f"[{time.time()-t0:.0f}s] Lam1..4 at 470 digits")
# cross-check vs integrals at ~150 digits (series identities F/H)
mp.mp.dps = 160
F = lambda a,b: mp.quad(lambda t: mp.log(t)*mp.log(1-b*t)/(1-a*t), [0,1])
H = lambda a,b: mp.quad(lambda t: mp.log(1-t)*mp.log(1-b*t)/(1-a*t), [0,1])
d1 = abs(F(x1,x2) - Lam['Lam1']/x1); d2 = abs(F(x2,x1) - Lam['Lam2']/x2)
d3 = abs(((-mp.log(1-x2))*mp.polylog(2,b1)/b1 + Lam['Lam3']/b1)/(1-x1) - H(x1,x2))
d4 = abs(((-mp.log(1-x1))*mp.polylog(2,b2)/b2 + Lam['Lam4']/b2)/(1-x2) - H(x2,x1))
print(f"cross-checks vs integrals (150dps): {mp.nstr(d1,3)} {mp.nstr(d2,3)} {mp.nstr(d3,3)} {mp.nstr(d4,3)}")
mp.mp.dps = 470
for k,v in Lam.items(): print(k, "=", mp.nstr(v, 60))

# ---------- 2. singles basis ----------
l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2 = mp.pi**2/6; z3 = mp.zeta(3)
logs = [(l2,'l2'),(l3,'l3'),(l5,'l5')]
def orbit(x):
    x = mp.mpf(x)
    fs = [lambda x:x, lambda x:1-x, lambda x:1/x, lambda x:x/(x-1), lambda x:1/(1-x), lambda x:(x-1)/x]
    out = set()
    for f in fs:
        try: out.add(mp.mpf(f(x)))
        except: pass
    return out
seeds = [x1, x2, mp.mpf(-13)/27, mp.mpf(-1)/5, mp.mpf(1)/6, mp.mpf(-7)/9, mp.mpf(-19)/81,
         mp.mpf(-1)/2, mp.mpf(-2)/3, mp.mpf(1)/2]
pts = set()
for s in seeds: pts |= orbit(s)
pts = sorted([p for p in pts if p < 1 and p != 0], key=lambda v: float(v))
print("Li points (%d):" % len(pts), [mp.nstr(p,6) for p in pts])
t0=time.time()
Li2v = {p: mp.polylog(2,p) for p in pts}
Li3v = {p: mp.polylog(3,p) for p in pts}
print(f"[{time.time()-t0:.0f}s] Li2/Li3 at {len(pts)} points")

basis = [(mp.mpf(1),'1')] + logs + [(z2,'z2')]
for p in pts: basis.append((Li2v[p], f'Li2({mp.nstr(p,8)})'))
for (a,na),(b,nb) in itertools.product(logs, repeat=2): basis.append((a*b, na+nb))
basis.append((z3,'z3'))
for p in pts: basis.append((Li3v[p], f'Li3({mp.nstr(p,8)})'))
for a,na in logs: basis.append((z2*a, 'z2'+na))
for p in pts:
    for a,na in logs: basis.append((Li2v[p]*a, f'Li2({mp.nstr(p,8)}){na}'))
for (a,na),(b,nb),(c,nc) in itertools.product(logs, repeat=3): basis.append((a*b*c, na+nb+nc))
print("basis size:", len(basis))
pickle.dump({'basis_names':[n for _,n in basis]}, open('basis_names.pkl','wb'))

# ---------- 3. LLL relation finder ----------
D = sp.Integer(10)**450
def find_relation(target, extra_note=""):
    """LLL on [basis, target]. Returns (verified_relations, exclusion_bound)."""
    vals = [v for v,_ in basis] + [target]
    n = len(vals)
    rows = []
    for i in range(n):
        r = [sp.Integer(0)]*n + [sp.Integer(int(mp.floor(vals[i]*D + mp.mpf('0.5'))))]
        r[i] = sp.Integer(1)
        rows.append(r)
    M = sp.Matrix(rows)
    L = M.lll(delta=sp.Rational(3,4))
    Lr = L.tolist()
    # candidate relations: rows with |last| < 1e250
    rels = []
    for r in Lr:
        if all(x==0 for x in r[:-1]): continue
        if abs(int(r[-1])) < 10**250:
            # verify residual at 470 digits
            res = mp.fsum([int(r[i])*vals[i] for i in range(n)])
            if abs(res) < mp.mpf(10)**(-430):
                rels.append((r, res))
    # exclusion bound via Gram-Schmidt of reduced basis
    # short vectors in lattice have first-n part = coefficients
    nonrel_firstn = min((sum(int(x)**2 for x in r[:-1]) for r in Lr if abs(int(r[-1])) >= 10**250), default=None)
    return rels, nonrel_firstn, n

results = {}
for nm, Lv in Lam.items():
    t0=time.time()
    rels, bound, n = find_relation(Lv, nm)
    dt=time.time()-t0
    tgt_rels = [(r,res) for r,res in rels if r[-1] != 0 or (r[n-1] != 0)]
    # target is at index n-1
    tgt_rels = [(r,res) for r,res in rels if r[n-1] != 0]
    if tgt_rels:
        r,res = tgt_rels[0]
        print(f"{nm}: REDUCES (relation found, residual {mp.nstr(abs(res),3)}, maxcoeff {max(abs(int(x)) for x in r[:-1])}) [{dt:.0f}s]")
        results[nm] = ('REDUCES', [int(x) for x in r])
    else:
        import math
        B = math.isqrt(bound) if bound else None
        print(f"{nm}: SURVIVES at 450dps over basis {n-1}+target; no relation with first-n norm^2 < {bound} (=> maxcoeff >~ {B}) [{dt:.0f}s]")
        results[nm] = ('SURVIVES', bound)
pickle.dump({'Lam': {k: mp.nstr(v, 460) for k,v in Lam.items()}, 'results': {k:(v[0], str(v[1])) for k,v in results.items()}},
            open('lam_lll_450.pkl','wb'))
with open("progress.log","a") as f:
    f.write(f"STEP 06: 470-digit constants; LLL tests: " + "; ".join(f"{k}:{v[0]}" for k,v in results.items()) + "\n")
print(f"total time {time.time()-t_start:.0f}s")
'''
with open("step_06.py","w") as f: f.write(code)
import subprocess, sys
r = subprocess.run([sys.executable, "step_06.py"], capture_output=True, text=True, cwd="research/08_theory/OCT7_generic_octagon", timeout=3000)
print(r.stdout)
if r.stderr: print("STDERR:", r.stderr[-800:])
