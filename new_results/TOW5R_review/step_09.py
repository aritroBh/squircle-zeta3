"""step_09 (TOW5R): sanity of the forced family + displayed closed forms + even-delta seed blindness.
 [1] F(dl,s) == tow3lib.e(dl,s) for dl <= 6, all s (exact).
 [2] F(dl, dl//2) == zigzag e_top(dl) for dl = 0..10 (exact).
 [3] displayed forms: F(9,4) = -16 i (d-6)/Pi_4; F(10,5) = -61/Pi_5; F(10,4) as in proof_depth10.md.
 [4] seed blindness (TOW4R caveat): perturb a_10 -> a_10 + eps (lower family fixed);
     d(F(10,6))/d eps == 0 exactly (so rho(10) rests on the R_5 derivation alone);
     and at odd delta = 9 the sensitivity is NONZERO (control).
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy')
import sympy as sp
import tow2lib as T2, tow3lib as T3, tow5lib as T5
d, m = T2.d, T2.m
i = sp.I
t0 = time.time()
ok1 = all(sp.cancel(sp.together(T5.forced(dl, s) - T3.e(dl, s))) == 0
          for dl in range(0, 7) for s in range(0, dl//2+1))
print('[1] F == tow3lib.e (dl<=6, all s):', ok1, flush=True)
ok2 = all(sp.cancel(sp.together(T5.forced(dl, dl//2) - T2.e_top(dl))) == 0 for dl in range(0, 11))
print('[2] F(top) == zigzag (dl<=10):', ok2, flush=True)
Pi4 = T5.Pi(4); Pi5 = T5.Pi(5)
f94 = -16*i*(d-6)/Pi4
f105 = -sp.Integer(61)/Pi5
f104 = (-16*d**3 + 288*d**2 + (-1777 + 8*i)*d - 16*i*m + 3554)/Pi5
ok3 = (sp.cancel(sp.together(T5.forced(9,4) - f94)) == 0 and
       sp.cancel(sp.together(T5.forced(10,5) - f105)) == 0 and
       sp.cancel(sp.together(T5.forced(10,4) - f104)) == 0)
print('[3] displayed F(9,4), F(10,5), F(10,4) match:', ok3, flush=True)

# [4] seed blindness: rebuild delta=10 closing with perturbed seed a_10 + eps
eps = sp.Symbol('eps')
import pickle
raw78 = pickle.load(open(os.path.dirname(os.path.abspath(__file__)) + '/tow5_copy/forced_F78.pkl','rb'))
Fz = {k: sp.sympify(v, locals={'d': d, 'm': m, 'I': sp.I}) for k, v in raw78.items()}
tau1 = T5.tau1f; tau2 = T5.tau2f
def make_FL(G):
    def F(dp, s):
        if s > dp//2: return sp.Integer(0)
        return G[(dp, s)]
    def Lam(delta, s):
        tot = sp.Integer(0)
        for dp in range(0, delta+1):
            fA = F(dp, s).subs({d: d+1, m: m+1}, simultaneous=True)
            fB = F(dp, s).subs({d: d+1, m: m}, simultaneous=True)
            tot += i*(2*d-1)/(2*d-dp-1) * (tau1(delta-dp)*fA - tau2(delta-dp)*fB)
        return tot
    return F, Lam
def build(delta, G):
    # build the unperturbed family at depth delta into G (if absent); return closing value
    F, Lam = make_FL(G)
    if (delta, 0) not in G:
        G[(delta, 0)] = sp.cancel(T5.arho(delta))
    val = sp.Integer(0)
    for s in range(0, delta//2+1):
        nxt = sp.together(Lam(delta, s) - F(delta, s))
        if s+1 <= delta//2:
            G[(delta, s+1)] = sp.cancel(nxt)
        else:
            val = nxt
    return val
def closing(delta, seed):
    G = dict(Fz)
    for dl in range(9, delta):
        build(dl, G)
    G[(delta, 0)] = seed
    F, Lam = make_FL(G)
    val = sp.Integer(0)
    for s in range(0, delta//2+1):
        nxt = sp.together(Lam(delta, s) - F(delta, s))
        if s+1 <= delta//2:
            G[(delta, s+1)] = sp.cancel(nxt)
        else:
            val = nxt
    return val
c10 = closing(10, sp.cancel(T5.arho(10)) + eps)
num, den = sp.fraction(sp.together(c10))
nume = sp.Poly(sp.expand(num), eps)
sens = sp.expand(nume.coeff_monomial(eps))
print('[4] delta=10: closing numerator as poly in eps: max eps-degree', sp.degree(num, eps),
      '; d(closing)/deps == 0 :', sens == 0, f'({round(time.time()-t0,1)}s)', flush=True)
c9 = closing(9, sp.cancel(T5.arho(9)) + eps)
num9 = sp.fraction(sp.together(c9))[0]
sens9 = sp.expand(sp.Poly(sp.expand(num9), eps).coeff_monomial(eps))
print('[4] delta=9 (control): sensitivity == 0 ?', sens9 == 0, '(expected NONZERO)', flush=True)
ok = ok1 and ok2 and ok3 and (sens == 0) and (sens9 != 0)
print('STEP_09 VERDICT:', 'PASS' if ok else 'FAIL', f'({round(time.time()-t0,1)}s)', flush=True)
assert ok
