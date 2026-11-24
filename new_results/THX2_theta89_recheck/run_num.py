from common import *
import mpmath as mp, itertools
mp.mp.dps = 75
Dtarget = pickle.load(open(T+'Dtilde_relations.pkl', 'rb'))
d89 = pickle.load(open(T+'I_reductions_w89.pkl', 'rb')); d9 = pickle.load(open(T+'I_reductions_w9.pkl', 'rb'))
basis = {8: d89['basis8'], 9: d89['basis9']}; reds = {8: d89['reds'], 9: d9}
def mpF(x): x = F(x); return mp.mpf(x.numerator)/x.denominator
# quadrature of atoms
Iq = {}
def I(r, k):
    if (r, k) not in Iq:
        Iq[(r, k)] = mp.quad(lambda s: s**r*mp.log(1+mp.e**(-s))**k, [0, 1, 4, 16, 64, 200])
    return Iq[(r, k)]
# own y-series evaluator for G(w;1), w in {0,-1}^n ending in -1 : t=y/(1-y), evaluate at y=1/2
NT = 330
half = mp.mpf(1)/2
def Gval(w):
    assert w[-1] == -1 and all(a in (0, -1) for a in w)
    f = [mp.mpf(1)] + [mp.mpf(0)]*NT
    for a in reversed(w):
        # g = int f * (1/(1-y) [+ 1/y if a==0]) dy
        run = mp.mpf(0); g = [mp.mpf(0)]*(NT+1)
        for n in range(NT):
            run += f[n]; term = run
            if a == 0: term += f[n+1]    # f/y contributes f_{n+1} y^n
            g[n+1] = term/(n+1)
        f = g
    return sum(f[k]*half**k for k in range(NT+1))
# (1) atom identity  I_{m,j} = j! m! G(0^{m+1}(-1)^j;1)  for all m+j+1<=9
worst = 0
for m in range(0, 8):
    for j in range(1, 9-m):
        if m+j+1 > 9: continue
        v = fact(j)*fact(m)*Gval((0,)*(m+1)+(-1,)*j)
        worst = max(worst, abs(v-I(m, j))/max(1, abs(v)))
print("atom identity I_{m,j}=j!m!G(0^{m+1}(-1)^j): worst abs error over all m+j+1<=9:", mp.nstr(worst, 3))
# (2) D~ vs zeta products
zv = lambda n: mp.zeta(n)
def zprod_val(ns):
    r = mp.mpf(1)
    for n in ns: r *= zv(n)
    return r
for W in range(5, 10):
    for m in range(W//2):
        D = sum(2*mpF(e_table(m, W-1-j))*I(W-1-j, j)/fact(j) for j in range(1, W-2*m))
        T_ = sum(mpF(c)*zprod_val(ns) for c, ns in zip(Dtarget[f"{W},{m}"], ZB[W]))
        print(f"D~_{{{W},{m}}} quadrature vs zeta products: |diff| = {mp.nstr(abs(D-T_), 3)}   (value {mp.nstr(D, 12)})")
# (3) basis symbol values
def li_half_vec(n):
    out = {}
    for c in itertools.product((0, -1), repeat=n-1):
        w = c + (-1,); out[w] = out.get(w, 0) + (-1)**sum(1 for a in c if a == -1)
    return out
for n in range(4, 9):
    v = sum(c*Gval(w) for w, c in li_half_vec(n).items())
    print(f"Li_{n}(1/2) via word combination vs polylog: {mp.nstr(abs(v-mp.polylog(n, half)), 3)}")
def zval(S):   # zeta of sum tuple with all word letters in {0,-1}: (-1)^d G(word)
    w = sum_to_word(S); return (-1)**len(S)*Gval(w)
X6 = mp.quad(lambda x: x**4*mp.log(1-mp.e**(-x))/(mp.e**x+1), [0, 1, 4, 16, 64, 200])/24
PV = {'log2': mp.log(2), 'X6': X6,
      'X7a': zval(((5,-1),(1,1),(1,1))), 'X7b': zval(((4,-1),(2,1),(1,1))),
      'X8a': zval(((6,-1),(1,1),(1,1))), 'X8b': zval(((5,-1),(2,1),(1,1))), 'X8c': zval(((5,-1),(1,1),(1,1),(1,1))), 'X8d': zval(((4,-1),(2,1),(1,1),(1,1))),
      'zbar(7,1,1)': zval(((7,-1),(1,1),(1,1))), 'zbar(6,2,1)': zval(((6,-1),(2,1),(1,1))), 'zbar(5,3,1)': zval(((5,-1),(3,1),(1,1))),
      'zbar(5,1,1,1,1)': zval(((5,-1),)+((1,1),)*4), 'zbar(4,2,1,1,1)': zval(((4,-1),(2,1),(1,1),(1,1),(1,1))),
      'X9fp': Gval((-1,-1,-1,0,0,0,-1,0,-1))}
def pval(nm):
    if nm in PV: return PV[nm]
    if nm[0] == 'z': return mp.zeta(int(nm[1:]))
    if nm.startswith('Li'): return mp.polylog(int(nm[2:]), half)
    raise KeyError(nm)
# the weight-6 reduction used X6 (independent integral value): check I_{3,2} reduction
r67 = pickle.load(open(T+'I_reductions_w2_7.pkl', 'rb'))
# (4) the weight 8/9 I-reductions
for W in (8, 9):
    names = basis[W]; vals = []
    for nm in names:
        p = mp.mpf(1)
        for f in nm.split('*'): p *= pval(f)
        vals.append(p)
    for m in range(W-1):
        j = W-1-m
        rhs = sum(mpF(c)*v for c, v in zip(reds[W][(m, j)], vals))
        lhs = I(m, j)
        print(f"I_{{{m},{j}}} (w{W}) quadrature vs reduction: |diff| = {mp.nstr(abs(lhs-rhs), 3)}  (|I|={mp.nstr(abs(lhs),6)})")
