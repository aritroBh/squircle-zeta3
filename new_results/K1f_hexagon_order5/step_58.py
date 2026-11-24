code = r'''
# step_26.py: numeric verification (>=60 digits) of the atom->HPL reduction table.
# Raw atom S((ms),a)(x) = sum_{n>=1} x^n prod_j (H_n^{(j)})^{ej}/n^a  vs HPL expression.
# HPL words evaluated as MPL single sums with exact inner finite sums.
import sympy as sp, mpmath as mp, pickle
mp.mp.dps = 80
WD = '<project>/research/08_theory/K1f_hexagon_order5'
tab = {k: sp.sympify(v) for k, v in pickle.load(open(f'{WD}/atom_table.pkl','rb')).items()}
H1, Li2, Li3, H011, Li4, H0011, H0111, Li5, H00011, H00101, H00111, H01011, H01111 = sp.symbols(
    'H1 Li2 Li3 H011 Li4 H0011 H0111 Li5 H00011 H00101 H00111 H01011 H01111')

# word -> MPL s-tuple
def word_to_mpl(w):
    s = []; run = 0
    for b in w:
        if b == 0: run += 1
        else: s.append(run+1); run = 0
    return tuple(s)

# numeric MPL  Li_{s1..sl}(x,1..1) = sum_{n1>...>nl>=1} x^n1 / prod ni^si
# evaluated as single sum over n1 with exact inner nested harmonic sums (recursion on depth).
def nested(ms, N):
    """sum_{N>j1>...>jl>=1} 1/prod ji^mi  computed exactly (mpf)."""
    if not ms: return mp.mpf(1)
    m0 = ms[0]
    tot = mp.mpf(0); inner = mp.mpf(0)
    # iterate j1 from 1 to N-1, accumulating inner for j2..
    # need inner(j1) = nested(ms[1:], j1)
    res = mp.mpf(0)
    cache = [mp.mpf(0)]*(N+1)
    if len(ms) > 1:
        # build incremental: nested(ms[1:], j) as j grows
        acc = mp.mpf(0)
        for j in range(1, N):
            acc += nested(ms[1:], j) / mp.mpf(j)**ms[0]
        res = acc
    else:
        for j in range(1, N):
            res += mp.mpf(1)/mp.mpf(j)**ms[0]
    return res

def mpl_num(stup, x, N=3000):
    s0 = stup[0]
    if len(stup) == 1:
        return mp.nsum(lambda n: x**n/mp.mpf(n)**s0, [1, mp.inf], tol=mp.mpf('1e-70'))
    return mp.nsum(lambda n: x**n * nested(stup[1:], n)/mp.mpf(n)**s0, [1, mp.inf], tol=mp.mpf('1e-60'))

# Lyndon atom word -> numeric
WORD2MPL = {w: word_to_mpl(w) for w in [(1,),(0,1),(0,0,1),(0,1,1),(0,0,0,1),(0,0,1,1),(0,1,1,1),
            (0,0,0,0,1),(0,0,0,1,1),(0,0,1,0,1),(0,0,1,1,1),(0,1,0,1,1),(0,1,1,1,1)]}
ATOMWORDS = {H1:(1,), Li2:(0,1), Li3:(0,0,1), H011:(0,1,1), Li4:(0,0,0,1), H0011:(0,0,1,1),
             H0111:(0,1,1,1), Li5:(0,0,0,0,1), H00011:(0,0,0,1,1), H00101:(0,0,1,0,1),
             H00111:(0,0,1,1,1), H01011:(0,1,0,1,1), H01111:(0,1,1,1,1)}

def atom_num(sym, x, N=3000):
    w = ATOMWORDS[sym]
    return mpl_num(WORD2MPL[w], x, N)

def raw_atom_num(a, ms, x):
    def term(n):
        v = mp.mpf(1)
        for j in set(ms):
            ej = ms.count(j)
            v *= (mp.fsum(mp.mpf(1)/mp.mpf(k)**j for k in range(1, n+1)))**ej
        return v * x**n / mp.mpf(n)**a
    return mp.nsum(term, [1, mp.inf], tol=mp.mpf('1e-60'))

import time
xs = [mp.mpf(1)/3, -mp.mpf(1)/2, mp.mpf(3)/4]
worst = mp.mpf(0); t0=time.time()
# first check the MPL sign convention on a few words
print("sign convention check Li_{2,1}(x,1) vs H_{0,1,1}(x)=S12-ish:")
x = mp.mpf(1)/3
v1 = mpl_num((2,1), x)
v2 = mp.nsum(lambda n: x**n * mp.fsum(mp.mpf(1)/k for k in range(1,n))/mp.mpf(n)**2, [1, mp.inf]) - mp.polylog(3, x)
print("  Li21 =", mp.nstr(v1,25), "  sum H_{n-1}/n^2 =", mp.nstr(v2,25), " diff", mp.nstr(abs(v1-v2),3))

results = []
for (a, ms), expr in sorted(tab.items()):
    msL = list(ms)
    for x in xs:
        lhs = raw_atom_num(a, msL, x)
        rhs = mp.mpf(0)
        for term in sp.Add.make_args(sp.expand(expr)):
            co, syms = term.as_coeff_mul()
            val = mp.mpf(1)
            for s_ in syms:
                if s_.is_Atom and s_.args:
                    base, ex = s_.args
                    val *= atom_num(base, x)**int(ex)
                else:
                    val *= atom_num(s_, x)
            rhs += mp.mpf(str(co)) * val
        d = abs(lhs-rhs)
        worst = max(worst, d)
        results.append(((a, ms), str(x), mp.nstr(d, 3)))
print("ALL atom reductions checked at 3 x-values, 80dps")
print("worst |lhs-rhs| over %d checks: %s" % (len(results), mp.nstr(worst, 3)))
for r in results:
    if mp.mpf(r[2]) > mp.mpf('1e-40'):
        print("  LARGE:", r)
print("time %.0f s" % (time.time()-t0))
'''
wd = '<project>/research/08_theory/K1f_hexagon_order5'
open(f'{wd}/step_26.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
subprocess.Popen([PY, f'{wd}/step_26.py'], cwd=wd, stdout=open(f'{wd}/step_26.out','w'), stderr=subprocess.STDOUT)
print("launched step_26 (numeric atom verification)")
