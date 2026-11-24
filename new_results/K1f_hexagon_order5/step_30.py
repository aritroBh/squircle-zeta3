
# step_30.py: FULL weight-5 non-local J-table build (210 entries) + numeric verification vs quad.
import sympy as sp, pickle, sys, time
WD = '<project>/research/08_theory/K1f_hexagon_order5'
sys.path.insert(0, WD)
import esum5
esum5.load_atom_table(WD)
x = esum5.x; Hf = esum5.Hf
m = sp.Symbol('m')
z2, z3, z5 = sp.symbols('z2 z3 z5')
jloc_raw = pickle.load(open(f'{WD}/jloc_sym.pkl','rb'))
Hs = sp.symbols('H1 H2 H3 H4 H5'); T2 = sp.symbols('T1 T2 T3 T4 T5'); DEN = sp.Symbol('DEN')
jloc = {k: sp.sympify(v) for k, v in jloc_raw.items()}

def jloc_sub(n, a, b):
    e = jloc[(a, b)]
    subs = {DEN: 1/(m+n+1)}
    for r in range(1, 6):
        subs[Hs[r-1]] = Hf(r, m+n+1)
        subs[T2[r-1]] = 1/(m+n+2)**r
    return sp.expand(e.subs(subs, simultaneous=True))

def esym(c, mv):
    H = [None] + [Hf(j, mv-1) for j in range(1, 5)]
    if c == 1: return sp.Integer(1)
    if c == 2: return H[1]
    if c == 3: return (H[1]**2 - H[2])/2
    if c == 4: return (H[1]**3 - 3*H[1]*H[2] + 2*H[3])/6
    if c == 5: return (H[1]**4 - 6*H[1]**2*H[2] + 3*H[2]**2 + 8*H[1]*H[3] - 6*H[4])/24
    raise ValueError(c)

def J(n, a, b, c):
    g = sp.expand(esym(c, m)/m * jloc_sub(n, a, b))
    return sp.expand(sp.factorial(c)*(-1)**c * esum5.esum(g, c, m))

t0 = time.time()
JTAB = {}
keys = [(n,a,b,c) for c in range(1,6) for a in range(6-c) for b in range(6-c-a) for n in range(6)]
print("total keys:", len(keys))
for i, (n,a,b,c) in enumerate(keys):
    JTAB[(n,a,b,c)] = J(n,a,b,c)
    if i % 30 == 0: print(f"{i}/{len(keys)} [{time.time()-t0:.0f}s]", flush=True)
pickle.dump({k: str(v) for k,v in JTAB.items()}, open(f'{WD}/jtab5.pkl','wb'))
print("saved jtab5.pkl, %.0f s" % (time.time()-t0))

# ---- numeric verification vs mpmath quad ----
import mpmath as mp, random
mp.mp.dps = 80
ATOMSYMS = sp.symbols('H1 Li2 Li3 H011 Li4 H0011 H0111 Li5 H00011 H00101 H00111 H01011 H01111')
def word_to_mpl(w):
    s = []; run = 0
    for b in w:
        if b == 0: run += 1
        else: s.append(run+1); run = 0
    return tuple(s)
ATOMWORDS = dict(zip(ATOMSYMS, [(1,),(0,1),(0,0,1),(0,1,1),(0,0,0,1),(0,0,1,1),(0,1,1,1),
            (0,0,0,0,1),(0,0,0,1,1),(0,0,1,0,1),(0,0,1,1,1),(0,1,0,1,1),(0,1,1,1,1)]))
NMAX = 2500
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
def evalJ(expr, xx):
    zv = {z2: mp.zeta(2), z3: mp.zeta(3), z5: mp.zeta(5)}
    tot = mp.mpf(0)
    for term in sp.Add.make_args(expr):
        co, fs = term.as_coeff_mul()
        val = mp.mpf(str(co))
        for f in fs:
            base, ex = (f.args if f.is_Pow else (f, 1))
            ex = int(ex)
            if base == x: val *= xx**ex
            elif base in zv: val *= zv[base]**ex
            elif base in ATOMWORDS: val *= atom_num(base, xx)**ex
            else: raise ValueError(f"unknown {f}")
        tot += val
    return tot
random.seed(2)
worst = mp.mpf(0); nchk = 0
kappas = [mp.mpf(1)/2, -mp.mpf(2)/3, mp.mpf(1)/3, -mp.mpf(1)/4]
tries = random.sample(keys, 40)
for (n,a,b,c) in tries:
    Je = JTAB[(n,a,b,c)]
    for kap in kappas:
        num = mp.quad(lambda t: t**n*mp.log(t)**a*mp.log(1-t)**b*mp.log(1+kap*t)**c, [0, mp.mpf('0.5'), 1])
        val = evalJ(Je, -kap)
        d = abs(num-val); worst = max(worst, d); nchk += 1
        if d > mp.mpf('1e-40'): print("  LARGE:", (n,a,b,c), kap, mp.nstr(d,3))
print("J-table checks: %d (40 keys x 4 kappas), worst |diff| = %s" % (nchk, mp.nstr(worst,3)))
with open(f'{WD}/progress.log','a') as f:
    f.write("step_30: full non-local J-table built (%d keys), saved jtab5.pkl; numeric verification vs quad: %d checks, worst |diff|=%s at 80dps.\n" % (len(keys), nchk, mp.nstr(worst,3)))
