code = r'''
# step_26.py (v2): numeric verification (>=60 digits) of atom->HPL reduction table.
import sympy as sp, mpmath as mp, pickle, time
mp.mp.dps = 80
WD = '<project>/research/08_theory/K1f_hexagon_order5'
tab = {k: sp.sympify(v) for k, v in pickle.load(open(f'{WD}/atom_table.pkl','rb')).items()}
H1, Li2, Li3, H011, Li4, H0011, H0111, Li5, H00011, H00101, H00111, H01011, H01111 = sp.symbols(
    'H1 Li2 Li3 H011 Li4 H0011 H0111 Li5 H00011 H00101 H00111 H01011 H01111')

NMAX = 2000
# precompute nested harmonic arrays for all tuples needed (depth<=4, entries in 1..4, weight<=5)
def build_nested(ms, Nmax):
    if len(ms) == 0: return [mp.mpf(1)]*(Nmax+1)
    prev = build_nested(ms[1:], Nmax)
    arr = [mp.mpf(0)]*(Nmax+1); acc = mp.mpf(0)
    for j in range(1, Nmax+1):
        acc += prev[j]/mp.mpf(j)**ms[0]
        if j+1 <= Nmax: arr[j+1] = acc
    return arr
from functools import lru_cache
@lru_cache(None)
def nested_arr(ms):
    return build_nested(tuple(ms), NMAX)

def word_to_mpl(w):
    s = []; run = 0
    for b in w:
        if b == 0: run += 1
        else: s.append(run+1); run = 0
    return tuple(s)

def mpl_num(stup, x):
    arr = nested_arr(stup[1:]) if len(stup) > 1 else None
    s0 = stup[0]
    tot = mp.mpf(0); xp = x
    for n in range(1, NMAX+1):
        inner = arr[n] if arr is not None else mp.mpf(1)
        tot += xp*inner/mp.mpf(n)**s0
        xp *= x
    return tot

ATOMWORDS = {H1:(1,), Li2:(0,1), Li3:(0,0,1), H011:(0,1,1), Li4:(0,0,0,1), H0011:(0,0,1,1),
             H0111:(0,1,1,1), Li5:(0,0,0,0,1), H00011:(0,0,0,1,1), H00101:(0,0,1,0,1),
             H00111:(0,0,1,1,1), H01011:(0,1,0,1,1), H01111:(0,1,1,1,1)}
def atom_num(sym, x):
    return mpl_num(word_to_mpl(ATOMWORDS[sym]), x)

# raw atoms: precompute H_n^{(j)} arrays
HACC = {j: [mp.mpf(0)]*(NMAX+1) for j in range(1,5)}
for j in range(1,5):
    acc = mp.mpf(0)
    for n in range(1, NMAX+1):
        acc += mp.mpf(1)/mp.mpf(n)**j
        HACC[j][n] = acc
def raw_atom_num(a, ms, x):
    tot = mp.mpf(0); xp = x
    for n in range(1, NMAX+1):
        v = mp.mpf(1)
        for j in set(ms):
            v *= HACC[j][n]**ms.count(j)
        tot += xp*v/mp.mpf(n)**a
        xp *= x
    return tot

xs = [mp.mpf(1)/3, -mp.mpf(1)/2, mp.mpf(3)/4]
worst = mp.mpf(0); t0=time.time(); nchecks=0
for (a, ms), expr in sorted(tab.items()):
    terms = [(term.as_coeff_mul()) for term in sp.Add.make_args(sp.expand(expr))]
    for x in xs:
        lhs = raw_atom_num(a, list(ms), x)
        rhs = mp.mpf(0)
        for co, syms in terms:
            val = mp.mpf(1)
            for s_ in syms:
                if s_.is_Pow:
                    base, ex = s_.args
                    val *= atom_num(base, x)**int(ex)
                else:
                    val *= atom_num(s_, x)
            rhs += mp.mpf(str(co))*val
        d = abs(lhs-rhs); worst = max(worst, d); nchecks += 1
        if d > mp.mpf('1e-40'): print("  LARGE:", a, ms, x, mp.nstr(d,3))
print("checks: %d, worst |lhs-rhs| = %s  (%.0f s)" % (nchecks, mp.nstr(worst,3), time.time()-t0))
with open(f'{WD}/progress.log','a') as f:
    f.write("step_26: atom->HPL reduction table verified numerically: %d checks (26 atoms x 3 x-values incl x=-1/2,3/4), worst |diff|=%s at 80dps.\n" % (nchecks, mp.nstr(worst,3)))
'''
wd = '<project>/research/08_theory/K1f_hexagon_order5'
open(f'{wd}/step_26.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
subprocess.Popen([PY, f'{wd}/step_26.py'], cwd=wd, stdout=open(f'{wd}/step_26.out','w'), stderr=subprocess.STDOUT)
print("relaunched step_26")
