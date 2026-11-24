
import mpmath as mp, sympy as sp, pickle, time, itertools, math
from fractions import Fraction
mp.mp.dps = 470
t_start = time.time()

# ---- constants at 470 ----
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
x1 = Fraction(-1,3); x2 = Fraction(-1,9); b1 = Fraction(1,4); b2 = Fraction(1,10)
def f2m(q): return mp.mpf(q.numerator)/q.denominator
Lam = {'Lam1': Li21(f2m(x1), f2m(x2/x1)), 'Lam2': Li21(f2m(x2), f2m(x1/x2)),
       'Lam3': Li21(f2m(b1), f2m(b2/b1)), 'Lam4': Li21(f2m(b2), f2m(b1/b2))}

# ---- E_arc at 470 digits via exact assembly ----
import mixed_engine_gen as E
Earcs = pickle.load(open('Earc_exact.pkl','rb'))['Earcs']
L2 = lambda v: mp.polylog(2, mp.mpf(v)); L3 = lambda v: mp.polylog(3, mp.mpf(v))
NUM = {
 E.l2: mp.log(2), E.l3: mp.log(3), E.l5: mp.log(5), E.z2: mp.pi**2/6, E.z3: mp.zeta(3),
 E.L2x1: L2(-mp.mpf(1)/3), E.L3x1: L3(-mp.mpf(1)/3), E.L2x2: L2(-mp.mpf(1)/9), E.L3x2: L3(-mp.mpf(1)/9),
 E.L2b1: L2(mp.mpf(1)/4), E.L3b1: L3(mp.mpf(1)/4), E.L2b2: L2(mp.mpf(1)/10), E.L3b2: L3(mp.mpf(1)/10),
 E.L2mh: L2(-mp.mpf(1)/2), E.L2m2t: L2(-mp.mpf(2)/3), E.L2m3h: L2(-mp.mpf(3)/2), E.L3m3h: L3(-mp.mpf(3)/2),
 E.L2m2: L2(-mp.mpf(2)), E.L3m2: L3(-mp.mpf(2)), E.L2_3f: L2(mp.mpf(3)/5), E.L3_3f: L3(mp.mpf(3)/5),
 E.L2_2t: L2(mp.mpf(2)/3), E.L3_2t: L3(mp.mpf(2)/3),
 E.Lam1: Lam['Lam1'], E.Lam2: Lam['Lam2'], E.Lam3: Lam['Lam3'], E.Lam4: Lam['Lam4'],
}
Earcs_num = mp.mpf(sp.N(Earcs.subs([(k, sp.Float(mp.nstr(v,470),470)) for k,v in NUM.items()]), 470))
print("E_arc at 470 digits:", mp.nstr(Earcs_num, 50))

# ---- independent basis at 470 ----
indep = pickle.load(open('basis_indep.pkl','rb'))['indep_names']
l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5)
z2, z3 = mp.pi**2/6, mp.zeta(3)
lm = {'l2':l2,'l3':l3,'l5':l5}
import re as _re
def bval(name):
    if name=='1': return mp.mpf(1)
    if name in lm: return lm[name]
    if name=='z2': return z2
    if name=='z3': return z3
    val = mp.mpf(1); rest = name
    if rest.startswith('z2'): val*=z2; rest=rest[2:]
    mm = _re.match(r'^Li([23])\((-?\d+/\d+|-?\d+)\)', rest)
    if mm:
        q = Fraction(mm.group(2)); val *= mp.polylog(int(mm.group(1)), f2m(q)); rest = rest[mm.end():]
    for ch in _re.findall(r'l[235]', rest): val *= lm[ch]
    return val
vals = [bval(nm) for nm in indep]
print("independent basis size:", len(vals))

# ---- LLL test ----
D = sp.Integer(10)**450
def test_target(target, label):
    valsL = vals + [target]
    n = len(valsL)
    rows = []
    for i in range(n):
        r = [sp.Integer(0)]*n + [sp.Integer(int(mp.floor(valsL[i]*D + mp.mpf('0.5'))))]
        r[i] = sp.Integer(1)
        rows.append(r)
    L = sp.Matrix(rows).lll(delta=sp.Rational(3,4)).tolist()
    # candidates: tiny last coord; verify at 470; then tolerance-stability at 700
    out = []
    for r in L:
        if all(x==0 for x in r[:-1]): continue
        if abs(int(r[-1])) < 10**240 and r[n-1] != 0:
            res = mp.fsum([int(r[i])*valsL[i] for i in range(n)])
            out.append((r, res))
    return out, n

results = {}
targets = list(Lam.items()) + [('E_arc', Earcs_num)]
for nm, tv in targets:
    t0=time.time()
    cands, n = test_target(tv, nm)
    good = [ (r,res) for r,res in cands if abs(res) < mp.mpf(10)**(-425) ]
    if good:
        r,res = good[0]
        print(f"{nm}: CANDIDATE relation, residual@470={mp.nstr(abs(res),3)}, maxcoeff {max(abs(int(x)) for x in r[:-1])} [{time.time()-t0:.0f}s] -> verify at 700")
        results[nm] = ('CANDIDATE', [int(x) for x in r])
    else:
        # exclusion: shortest candidate residual scale
        print(f"{nm}: SURVIVES over independent basis {n-1} at 450 digits (no target relation with residual < 1e-425) [{time.time()-t0:.0f}s]")
        results[nm] = ('SURVIVES', None)
pickle.dump({'results': results, 'Earc470': mp.nstr(Earcs_num,460)}, open('test_indep_450.pkl','wb'))
with open("progress.log","a") as f:
    f.write("STEP 08 (indep basis 90, 450d): " + "; ".join(f"{k}:{v[0]}" for k,v in results.items()) + "\n")
print(f"total {time.time()-t_start:.0f}s")
