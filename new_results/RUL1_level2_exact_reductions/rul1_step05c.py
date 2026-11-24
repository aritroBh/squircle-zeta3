
import sys, time, pickle
sys.dont_write_bytecode = True
import sympy as sp
from fractions import Fraction as Fr
import rul1, rul1_reduce

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t0 = time.time()
log("rul1_step05c: conv-B fold (unit normals, log3 tracked) vs E5e's X")

D = rul1_reduce.load()
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
coef_acc = pickle.load(open(f'{K1G}/coef_acc_raw.pkl','rb'))
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
l2, l3, pi2, z3, z5, Li4h, Li5h, sq5 = sp.symbols('l2 l3 pi2 z3 z5 Li4h Li5h sq5')

ATOMWORDS = {'H1':(1,), 'Li2':(0,1), 'Li3':(0,0,1), 'H011':(0,1,1), 'Li4':(0,0,0,1),
             'H0011':(0,0,1,1), 'H0111':(0,1,1,1), 'Li5':(0,0,0,0,1), 'H00011':(0,0,0,1,1),
             'H00101':(0,0,1,0,1), 'H00111':(0,0,1,1,1), 'H01011':(0,1,0,1,1), 'H01111':(0,1,1,1,1)}
def atom_combo(nm):
    hw = ATOMWORDS[nm]
    gw = tuple(-1 if a == 1 else 0 for a in hw)
    return {gw: Fr((-1)**sum(1 for a in hw if a == 1))}
ATOMW = {nm: len(w) for nm, w in ATOMWORDS.items()}
ATOMC = {nm: atom_combo(nm) for nm in ATOMWORDS}
ZETAC = {0: {(0,1): Fr(-1)}, 1: {(0,0,1): Fr(-1)}, 2: {(0,0,0,0,1): Fr(-1)}}
ZW = {0: 2, 1: 3, 2: 5}
from functools import lru_cache
@lru_cache(None)
def shpow(w, n):
    out = {(): Fr(1)} if n == 0 else {w: Fr(1)}
    for _ in range(n-1):
        out = rul1.sh_mul(out, {w: Fr(1)})
    return out
@lru_cache(None)
def reduce_cached(w, combo_key):
    if w == 0:
        d = dict(combo_key); assert set(d) <= {()}
        return (('1', d.get((), Fr(0))),)
    return tuple(sorted(rul1_reduce.reduce_combo(w, dict(combo_key)).items()))
def namedsym_expr(coords):
    expr = 0
    for nm, x in coords:
        if nm == '1': expr += x
        elif nm == 'pi4': expr += x*pi2**2
        elif nm == 'log2': expr += x*l2
        elif nm.startswith('log2^'): expr += x*l2**int(nm[5:])
        elif nm == 'pi2': expr += x*pi2
        elif nm == 'z3*pi2': expr += x*z3*pi2
        elif nm == 'z3': expr += x*z3
        elif nm == 'z5': expr += x*z5
        elif nm == 'Li4h': expr += x*Li4h
        elif nm == 'Li5h': expr += x*Li5h
        elif nm == 'pi2*log2^2': expr += x*pi2*l2**2
        elif nm == 'pi2*log2^3': expr += x*pi2*l2**3
        elif nm == 'z3*log2': expr += x*z3*l2
        elif nm == 'z3*log2^2': expr += x*z3*l2**2
        elif nm == 'pi2*log2': expr += x*pi2*l2
        elif nm == 'pi4*log2': expr += x*pi2**2*l2
        elif nm == 'Li4h*log2': expr += x*Li4h*l2
        else: raise ValueError(nm)
    return expr

def fold(rd):
    total = sp.Integer(0)
    nred = 0
    for (at, xe, ze), co in coef_acc.items():
        cov = sp.expand(co.subs(rd))
        if cov == 0: continue
        base = {(): Fr(1)}
        zw = 0
        for nm, e in at:
            for _ in range(e):
                base = rul1.sh_mul(base, ATOMC[nm])
            zw += ATOMW[nm]*e
        for i, e in enumerate(ze):
            for _ in range(e):
                base = rul1.sh_mul(base, ZETAC[i])
            zw += ZW[i]*e
        poly = sp.Poly(cov, l2, l3) if (cov.has(l2) or cov.has(l3)) else None
        terms = poly.terms() if poly else [((0,0), cov)]
        for (m2, m3), cm in terms:
            combo = dict(base)
            if m2:
                combo = rul1.sh_mul(combo, shpow((-1,), m2))
            wtot = zw + m2
            red = reduce_cached(wtot, tuple(sorted(combo.items())))
            nred += 1
            total += cm*(-1)**xe*namedsym_expr(red)*(l3**m3 if m3 else 1)
    return sp.expand(120*total), nred

# conv-B (unit-normal) logs: log(A_j/|m_j|) = (log2, -log2/2, 0); no log3 anywhere
rdB = {A0:sp.Rational(1,2),A1:sp.Rational(3,4),B:sp.Rational(1,4),kap:sp.Integer(1),
      L2:sp.Integer(5),q0:sp.Integer(1),V02:sp.Integer(1),lA0:l2,lA1:-l2/2,lB:0}
valB, nred = fold(rdB)
valB = valB.subs({sp.sqrt(5): sq5})
log("conv-B fold done (%d reductions, %.1fs)" % (nred, time.time()-t0))
# check log3 cancels
import collections
c = collections.Counter()
for t in sp.Add.make_args(sp.expand(valB)):
    p3 = sp.degree(t, l3) if t.has(l3) else 0
    c[p3] += 1
print("terms by l3-power:", dict(c))
l3part = sp.expand(valB).coeff(l3, 1)
assert not any(k for k in c if k > 0), "unexpected log3"
log("conv-B (unit-normal) value is log3-free: unit-normal logs are (log2, -log2/2, 0)")

# E5e closed form for X = P_arc/sqrt5
X_e5e = (-153*l2 - sp.Rational(32,5)*pi2 - sp.Rational(26817,80)*l2**2 + sp.Rational(3267,50)*z3
         + sp.Rational(291,20)*pi2*l2 - sp.Rational(386517,3200)*l2**3 - sp.Rational(6489,8000)*pi2**2
         - sp.Rational(20679,400)*z3*l2 - sp.Rational(18823,3200)*pi2*l2**2 - sp.Rational(382819,10240)*l2**4
         - sp.Rational(1242,25)*Li4h - sp.Rational(20187,2000)*z5 - sp.Rational(858,125)*pi2*z3/6
         + sp.Rational(6463,40000)*pi2**2*l2 + sp.Rational(57,250)*z3*l2**2 + sp.Rational(31267,24000)*pi2*l2**3
         + sp.Rational(992499,2560000)*l2**5 + sp.Rational(396,125)*Li5h + sp.Rational(918,125)*Li4h*l2)
diff = sp.simplify(sp.expand(valB - sq5*X_e5e))
print("valB - sqrt5*X_e5e =", diff)
assert diff == 0
log("TARGET 2c PASS: conv-B fold == sqrt5 * E5e's closed form X, EXACTLY (all 20 constants)")

# numeric: vs E5e's 336-digit value
import mpmath as mp
mp.mp.dps = 120
NV = rul1num_ = None
import rul1num
NV = rul1num.named_values(120)
val_num = valB.subs({l2: NV['log2'], pi2: NV['pi2'], z3: NV['z3'], z5: NV['z5'],
                     Li4h: NV['Li4h'], Li5h: NV['Li5h'], sq5: mp.sqrt(5)})
val_num = mp.mpf(str(val_num))
Xstr = "-381.140137469830729527176573357265848015393511891093057029373066042680677255274504566066309782740023456"
Xe5e = mp.mpf(Xstr)*mp.sqrt(5)
print("fold value      =", mp.nstr(val_num, 40))
print("E5e (336 digit) =", mp.nstr(Xe5e, 40))
print("|diff| =", mp.nstr(abs(val_num - Xe5e), 3))
pickle.dump(valB, open('rul1_fold_edge1_convB.pkl','wb'))
log("saved rul1_fold_edge1_convB.pkl, %.1fs" % (time.time()-t0))
