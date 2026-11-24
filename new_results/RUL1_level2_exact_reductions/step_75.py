import pathlib
code = r'''
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
log("rul1_step05b: exact fold of K1g's 307-monomial assembly at edge 1 (r=2, x=-1)")

D = rul1_reduce.load()
K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
coef_acc = pickle.load(open(f'{K1G}/coef_acc_raw.pkl','rb'))

A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
l2, l3, pi2, z3, z5, Li4h, Li5h, sq5 = sp.symbols('l2 l3 pi2 z3 z5 Li4h Li5h sq5')

ATOMWORDS = {'H1':(1,), 'Li2':(0,1), 'Li3':(0,0,1), 'H011':(0,1,1), 'Li4':(0,0,0,1),
             'H0011':(0,0,1,1), 'H0111':(0,1,1,1), 'Li5':(0,0,0,0,1), 'H00011':(0,0,0,1,1),
             'H00101':(0,0,1,0,1), 'H00111':(0,0,1,1,1), 'H01011':(0,1,0,1,1), 'H01111':(0,1,1,1,1)}
ATOMW = {nm: sum(1 for a in w if a == 1) and len(w) or len(w) for nm, w in ATOMWORDS.items()}
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

def combo_pow(combo, n):
    out = {(): Fr(1)}
    key = tuple(sorted(combo.items()))
    for _ in range(n):
        out = rul1.sh_mul(out, combo)
    return out

@lru_cache(None)
def reduce_cached(w, combo_key):
    return tuple(sorted(rul1_reduce.reduce_combo(w, dict(combo_key)).items()))

def namedsym_expr(coords):
    expr = 0
    for nm, x in coords:
        if nm == 'pi4': expr += x*pi2**2
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
        elif nm == 'Li4h*log2': expr += x*Li4h*l2
        else: raise ValueError(nm)
    return expr

def fold(rd, logsymbols):
    """rd: substitution dict for A0..lB (logs -> exprs in logsymbols).
       Returns total E_arc value as sympy expr in named constants + logsymbols beyond l2."""
    total = sp.Integer(0)
    nred = 0
    for (at, xe, ze), co in coef_acc.items():
        cov = sp.expand(co.subs(rd))
        if cov == 0: continue
        # monomial word-combo for atoms and zetas
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
        # coefficient as polynomial in the log symbols
        logsym = [l2] + ([l3] if l3 in cov.free_symbols else [])
        poly = sp.Poly(cov, *logsym) if cov.has(l2) or (l3 in cov.free_symbols) else None
        terms = poly.terms() if poly else [(tuple([0]*len(logsym)), cov)]
        for powers, cm in terms:
            m2 = powers[0]; m3 = powers[1] if len(logsym) > 1 else 0
            combo = dict(base)
            if m2:
                combo = rul1.sh_mul(combo, shpow((-1,), m2))
            wtot = zw + m2
            red = reduce_cached(wtot, tuple(sorted(combo.items())))
            nred += 1
            total += cm*(-1)**xe*namedsym_expr(red)*(l3**m3 if m3 else 1)
    return sp.expand(120*total), nred

rd1 = {A0:sp.Rational(1,2),A1:sp.Rational(3,4),B:sp.Rational(1,4),kap:sp.Integer(1),
      L2:sp.Integer(5),q0:sp.Integer(1),V02:sp.Integer(1),lA0:-l2,lA1:-2*l2,lB:-2*l2}
val, nred = fold(rd1, [l2])
val = val.subs({sp.sqrt(5): sq5})
log("polar fold done (%d reductions, %.1fs)" % (nred, time.time()-t0))
val = sp.expand(val)
cz5 = sp.expand(val).coeff(z5, 1)
cpz3 = sp.expand(val).coeff(pi2, 1).coeff(z3, 1) if val.has(pi2) else 0
# safer: extract from terms
terms = sp.Add.make_args(val)
def coeff_of(symexpr):
    s = 0
    for t in terms:
        cf = t/symexpr
        if not cf.has(z3) and not cf.has(z5) and not cf.has(Li4h) and not cf.has(Li5h) and not cf.has(pi2) and not cf.has(l2):
            s += t
    return sp.simplify(s/symexpr)
cz5 = coeff_of(z5)
cpi2z3 = coeff_of(pi2*z3)
print("z5 coefficient of E_arc      =", cz5)
print("pi2*z3 coefficient of E_arc  =", cpi2z3)
print("expected: -20187/2000*sqrt5  =", -sp.Rational(20187,2000)*sq5)
print("expected: z2z3 = pi2 z3/6: 6*(pi2z3 coef) =", sp.simplify(6*cpi2z3), " vs -858/125 sqrt5 =", -sp.Rational(858,125)*sq5)
ok1 = sp.simplify(cz5 + sp.Rational(20187,2000)*sq5) == 0
ok2 = sp.simplify(6*cpi2z3 + sp.Rational(858,125)*sq5) == 0
print("z5 match:", ok1, "  z2z3 match:", ok2)
assert ok1 and ok2
# non-local parts (subtract local edge law, X = P/sqrt5 normalization)
nl5 = sp.simplify(cz5 - (-sp.Rational(2079,250))*sq5)
nl23 = sp.simplify(6*cpi2z3 - (-sp.Rational(477,250))*sq5)
print("non-local z5   = total - (-2079/250 sqrt5) =", nl5, "  expect -711/400 sqrt5 =", -sp.Rational(711,400)*sq5)
print("non-local z2z3 = total - (-477/250 sqrt5) =", nl23, "  expect -1239/250 sqrt5 =", -sp.Rational(1239,250)*sq5)
assert sp.simplify(nl5 + sp.Rational(711,400)*sq5) == 0
assert sp.simplify(nl23 + sp.Rational(1239,250)*sq5) == 0
log("TARGET 2b PASS: exact fold gives z5 = -20187/2000 sqrt5, z2z3 = -858/125 sqrt5; far-kink parts -711/400 sqrt5, -1239/250 sqrt5")
import pickle
pickle.dump(val, open('rul1_fold_edge1_polar.pkl','wb'))
log("saved rul1_fold_edge1_polar.pkl, %.1fs" % (time.time()-t0))
'''
pathlib.Path('rul1_step05b.py').write_text(code)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice',PY,'rul1_step05b.py'], capture_output=True, text=True, timeout=3000)
print(r.stdout[-3000:]); print(r.stderr[-1000:])
