
# step_28.py: self-checks for esum5 engine at 80 dps against direct numeric sums.
import sympy as sp, mpmath as mp, pickle, sys, time
mp.mp.dps = 80
WD = '<project>/research/08_theory/K1f_hexagon_order5'
sys.path.insert(0, WD)
import esum5
esum5.load_atom_table(WD)
x = esum5.x
H1,H2,H3,H4 = esum5.HSYMS[1], esum5.HSYMS[2], esum5.HSYMS[3], esum5.HSYMS[4]
Hf = esum5.Hf
k = sp.Symbol('k')

ATOMSYMS = sp.symbols('H1 Li2 Li3 H011 Li4 H0011 H0111 Li5 H00011 H00101 H00111 H01011 H01111')

def word_to_mpl(w):
    s = []; run = 0
    for b in w:
        if b == 0: run += 1
        else: s.append(run+1); run = 0
    return tuple(s)
ATOMWORDS = dict(zip(ATOMSYMS, [(1,),(0,1),(0,0,1),(0,1,1),(0,0,0,1),(0,0,1,1),(0,1,1,1),
            (0,0,0,0,1),(0,0,0,1,1),(0,0,1,0,1),(0,0,1,1,1),(0,1,0,1,1),(0,1,1,1,1)]))

NMAX = 2000
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
def nested_arr(ms): return build_nested(tuple(ms), NMAX)
def mpl_num(stup, xx):
    s0 = stup[0]
    arr = nested_arr(stup[1:]) if len(stup) > 1 else None
    tot = mp.mpf(0); xp = xx
    for n in range(1, NMAX+1):
        inner = arr[n] if arr is not None else mp.mpf(1)
        tot += xp*inner/mp.mpf(n)**s0
        xp *= xx
    return tot
def atom_num(sym, xx): return mpl_num(word_to_mpl(ATOMWORDS[sym]), xx)

def eval_esum_expr(expr, xx):
    zv = {sp.Symbol('z2'): mp.zeta(2), sp.Symbol('z3'): mp.zeta(3), sp.Symbol('z5'): mp.zeta(5)}
    tot = mp.mpf(0)
    for term in sp.Add.make_args(sp.expand(expr)):
        co, fs = term.as_coeff_mul()
        val = mp.mpf(str(co))
        for f in fs:
            if f.is_Pow:
                base, ex = f.args; ex = int(ex)
            else:
                base, ex = f, 1
            if base == x:
                val *= xx**ex
            elif base in zv:
                val *= zv[base]**ex
            elif base in ATOMWORDS:
                val *= atom_num(base, xx)**ex
            elif base.funcname if hasattr(base,'funcname') else False:
                raise ValueError(f"unknown factor {f}")
            else:
                raise ValueError(f"unknown factor {f} in {term}")
        tot += val
    return tot

def num_sum(gfun, K, xx, N=4500):
    tot = mp.mpf(0); xp = xx**K
    for n in range(K, K+N):
        tot += xp*mp.mpf(gfun(n))
        xp *= xx
    return tot

NMAXN = 5000
Hnum = {j: [mp.mpf(0)]*(NMAXN+1) for j in range(1,5)}
for j in range(1,5):
    acc = mp.mpf(0)
    for n in range(1, NMAXN+1):
        acc += mp.mpf(1)/mp.mpf(n)**j
        Hnum[j][n] = acc
H = lambda n,j: Hnum[j][n]
tests = [
    ("sum x^k Hk/k^2", Hf(1,k)/k**2, 1, lambda n: H(n,1)/n**2),
    ("sum x^k Hk/(k+2)^2", Hf(1,k)/(k+2)**2, 1, lambda n: H(n,1)/(n+2)**2),
    ("sum x^k Hk^2 H2k/(k+1)", Hf(1,k)**2*Hf(2,k)/(k+1), 1, lambda n: H(n,1)**2*H(n,2)/(n+1)),
    ("sum x^k Hk/(k(k+1))", Hf(1,k)/(k*(k+1)), 1, lambda n: H(n,1)/(n*(n+1))),
    ("sum x^k Hk^3/(k+3)", Hf(1,k)**3/(k+3), 1, lambda n: H(n,1)**3/(n+3)),
    ("sum x^k Hk H3k/(k+2)", Hf(1,k)*Hf(3,k)/(k+2), 1, lambda n: H(n,1)*H(n,3)/(n+2)),
    ("sum x^k Hk^2/(k+1)^3", Hf(1,k)**2/(k+1)**3, 1, lambda n: H(n,1)**2/(n+1)**3),
    ("sum_{k>=3} x^k Hk^2/(k-1)", Hf(1,k)**2/(k-1), 3, lambda n: H(n,1)**2/(n-1)),
    ("sum x^k H4k/k", Hf(4,k)/k, 1, lambda n: H(n,4)/n),
]
xs = [mp.mpf(1)/3, -mp.mpf(1)/2, mp.mpf(3)/5]
worst = mp.mpf(0)
t0 = time.time()
for name, g, K, gn in tests:
    red = esum5.esum(g, K, k)
    for xx in xs:
        lhs = num_sum(gn, K, xx)
        rhs = eval_esum_expr(red, xx)
        d = abs(lhs-rhs); worst = max(worst, d)
        print(f"{name:32s} x={mp.nstr(xx,4)}: num={mp.nstr(lhs,12)}  esum={mp.nstr(rhs,12)}  |d|={mp.nstr(d,2)}")
print("worst over all:", mp.nstr(worst, 3), "  time %.0f s" % (time.time()-t0))
