import mpmath as mp, pickle, time
from fractions import Fraction
t0=time.time()
# ============ PART A: verify the 450-digit candidates at 700 digits ============
mp.mp.dps = 700
def Li21(a, c, digits=690):
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
Lam700 = {'Lam1': Li21(f2m(x1), f2m(x2/x1)), 'Lam2': Li21(f2m(x2), f2m(x1/x2)),
          'Lam3': Li21(f2m(b1), f2m(b2/b1)), 'Lam4': Li21(f2m(b2), f2m(b1/b2))}
l2, l3, l5 = mp.log(2), mp.log(3), mp.log(5); z2, z3 = mp.pi**2/6, mp.zeta(3)
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
indep = pickle.load(open('basis_indep.pkl','rb'))['indep_names']
vals700 = [bval(nm) for nm in indep]
print(f"[{time.time()-t0:.0f}s] 700-digit constants done")
d = pickle.load(open('test_indep_450.pkl','rb'))
# E_arc at 700 via exact coefficients: E_arc = sum over named constants; use stored exact assembly
import sympy as sp, mixed_engine_gen as E
Earcs = pickle.load(open('Earc_exact.pkl','rb'))['Earcs']
L2f = lambda v: mp.polylog(2, mp.mpf(v)); L3f = lambda v: mp.polylog(3, mp.mpf(v))
NUM = {E.l2:l2,E.l3:l3,E.l5:l5,E.z2:z2,E.z3:z3,
 E.L2x1:L2f(-mp.mpf(1)/3),E.L3x1:L3f(-mp.mpf(1)/3),E.L2x2:L2f(-mp.mpf(1)/9),E.L3x2:L3f(-mp.mpf(1)/9),
 E.L2b1:L2f(mp.mpf(1)/4),E.L3b1:L3f(mp.mpf(1)/4),E.L2b2:L2f(mp.mpf(1)/10),E.L3b2:L3f(mp.mpf(1)/10),
 E.L2mh:L2f(-mp.mpf(1)/2),E.L2m2t:L2f(-mp.mpf(2)/3),E.L2m3h:L2f(-mp.mpf(3)/2),E.L3m3h:L3f(-mp.mpf(3)/2),
 E.L2m2:L2f(-mp.mpf(2)),E.L3m2:L3f(-mp.mpf(2)),E.L2_3f:L2f(mp.mpf(3)/5),E.L3_3f:L3f(mp.mpf(3)/5),
 E.L2_2t:L2f(mp.mpf(2)/3),E.L3_2t:L3f(mp.mpf(2)/3),
 E.Lam1:Lam700['Lam1'],E.Lam2:Lam700['Lam2'],E.Lam3:Lam700['Lam3'],E.Lam4:Lam700['Lam4']}
Earc700 = mp.mpf(sp.N(Earcs.subs([(k, sp.Float(mp.nstr(v,700),700)) for k,v in NUM.items()]), 700))
tgt700 = dict(Lam700); tgt700['E_arc'] = Earc700
print(f"[{time.time()-t0:.0f}s] E_arc 700 done")
verdicts = {}
for nm, v in d['results'].items():
    if v[0] != 'CANDIDATE': verdicts[nm] = ('SURVIVES', None); continue
    r = v[1]  # full row: 90 basis coeffs + target coeff + value column
    res = mp.fsum([mp.mpf(r[i])*vals700[i] for i in range(90)]) + mp.mpf(r[90])*tgt700[nm]
    verdict = "GENUINE" if abs(res) < mp.mpf(10)**(-650) else "SPURIOUS"
    print(f"{nm}: candidate (maxcoeff {max(abs(x) for x in r[:90])}) residual@700 = {mp.nstr(abs(res),3)} -> {verdict}")
    verdicts[nm] = (verdict, mp.nstr(abs(res),3))
with open("progress.log","a") as f:
    f.write("STEP 08c @700: " + "; ".join(f"{k}:{v[0]}" for k,v in verdicts.items()) + "\n")
print(f"[{time.time()-t0:.0f}s] PART A done")
pickle.dump(verdicts, open('verdicts700.pkl','wb'))
