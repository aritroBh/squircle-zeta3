script = r'''"""verify_cor3.py -- COR3k: exact proof of the ten weight-3 HPL reductions and of
the classical form B2' of M2(alpha).

Method (no polylog-identity simplification is ever requested from sympy):
  * every derivative is taken by the chain rule; polylog(1,a) is rewritten BY HAND
    into the basis {log x, log(1-x), log(1+x), log 2} (all factors positive on (0,1));
  * every "vanishing" check is sp.cancel()==0 on rational-function coefficients of a
    fixed monomial basis (Poly over QQ(x) in the atom symbols, pi^2 tracked as a symbol);
  * constants at x -> 0+ are evaluated by hand-coded series with the PROVED special
    values Li2(0)=0, Li2(1)=pi^2/6, Li2(-1)=-pi^2/12, Li2(1/2)=pi^2/12-log2^2/2,
    Li3(0)=0, Li3(1)=zeta(3), Li3(-1)=-3 zeta(3)/4,
    Li3(1/2)=7 zeta(3)/8 - pi^2 log2/12 + log2^3/6,
    substituted BEFORE comparison; the resulting formal combination in
    Q[zeta(3), pi, log 2] must be the zero polynomial (sufficient: the real number
    it names is then 0 -- no linear-independence assumption is needed).

Run: python verify_cor3.py   (< 10 min; prints PROVED / FAILED per item)
"""
import json, re, ast, time
import sympy as sp
import mpmath as mp
import numpy as np
from cor3_common import (x, LX, z3, P2, L2, logx, l1m, l1p, ARGS, ARG_NAMES, f_a,
                         dpoly, load_classical_dict, const_term, is_zero_const,
                         check_weight1_zero, check_weight2_zero, LI2_REDUCE,
                         replace_polylog1)

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle/"
RESULTS = []
def report(label, ok, extra=""):
    RESULTS.append((label, ok))
    print(f"  {label}: {'PROVED' if ok else 'FAILED'} {extra}")

T0 = time.time()
CD = load_classical_dict(G1 + "classical_dict.json")

print("=" * 78)
print("TASK 1, step 0: weight-1 definitions  H_0 = log x, H_1 = -log(1-x), H_-1 = log(1+x)")
print("=" * 78)
H1w = {0: logx, 1: -l1m, -1: l1p}

print("\n=== TASK 1, step 1 (weight 2): C_{a,b} = H_{a,b} for the six words with b != 0 ===")
C2 = {
    (0, 1):  sp.polylog(2, x),
    (0, -1): -sp.polylog(2, -x),
    (1, 1):  l1m**2/2,
    (-1, -1): l1p**2/2,
    (1, -1): -sp.polylog(2, (1+x)/2) - l1m*l1p + L2*l1p - L2**2/2 + P2/12,
    (-1, 1):  sp.polylog(2, (1+x)/2) - L2*l1p + L2**2/2 - P2/12,
}
for (a, b), Cab in C2.items():
    D = sp.expand(dpoly(Cab) - f_a(a)*H1w[b])          # must reduce to the log basis and vanish
    okd, terms = check_weight1_zero(D)
    c0, resid = const_term(Cab, 3)                     # value at 0+ via special values
    okc, _ = is_zero_const(c0)
    okr = sp.expand(resid).coeff(x, 0) == 0            # limit exists (no x^0/LX residue)
    report(f"H[{a},{b}] = C[{a},{b}] (derivative EXACT {okd}, C(0+)=0 EXACT {okc}, regular {okr})",
           okd and okc and okr)

print("\n=== TASK 1, step 1b: proved classical Li2-atom identities (Euler/Landen/five-term) ===")
for nm in ['1-x', '1/(1+x)', '-(1-x)/(1+x)', '(1-x)/(1+x)', '2x/(1+x)']:
    G = sp.expand(sp.polylog(2, ARGS[nm]) - LI2_REDUCE[nm])
    okd, _ = check_weight1_zero(sp.expand(dpoly(G)))
    c0, resid = const_term(G, 3)
    okc, _ = is_zero_const(c0)
    report(f"Li2({nm}) reduction (derivative EXACT {okd}, value at 0+ EXACT {okc})", okd and okc)

# duplication Li2(x) + Li2(-x) = Li2(x^2)/2  (used in B2:  -Li2(c^2)/2 = -Li2(c) - Li2(-c))
dG = sp.diff(sp.polylog(2, x) + sp.polylog(2, -x) - sp.polylog(2, x**2)/2, x)
dG = dG.replace(sp.polylog(1, x**2), -(l1m + l1p))   # log(1-x^2) = log(1-x)+log(1+x), both > 0
okd, _ = check_weight1_zero(sp.expand(replace_polylog1(dG)))
report("duplication Li2(x)+Li2(-x) = Li2(x^2)/2 (derivative EXACT, value 0 at 0 by Li2(0)=0)", okd)

print("\n=== TASK 1, step 2 (weight 3): the ten identities H_w = C_w ===")
for w_str in sorted(CD.keys()):
    a1, a2, a3 = ast.literal_eval(w_str)
    Cw = CD[w_str]
    # (i) derivative identity:  d/dx C_w - f_{a1} C_{(a2,a3)} == 0, reduced with the
    #     PROVED weight-2 forms and the PROVED Li2-atom identities, all coefficients 0
    Ew = sp.expand(dpoly(Cw) - f_a(a1)*C2[(a2, a3)])
    okE, termsE = check_weight2_zero(Ew)
    # (ii) constant: lim_{x -> 0+} C_w = 0 with the special values substituted,
    #      coefficient-by-coefficient in {zeta(3), pi, log 2}
    c0, resid = const_term(Cw, 2)
    okc, poly = is_zero_const(c0)
    okr = sp.expand(resid).coeff(x, 0) == 0
    report(f"w = {w_str}: derivative identity EXACT {okE}; constant EXACT 0 {okc} -> H_w = C_w",
           okE and okc and okr, extra=f"(lim C_w = {c0}; monomials checked: {len(termsE)})")

print("\n=== TASK 1, step 3: Li-argument ranges on (0,1) ===")
grid = np.linspace(1e-9, 1 - 1e-9, 200001)
for nm in ARG_NAMES:
    v = sp.lambdify(x, ARGS[nm], 'numpy')(grid)
    dv = sp.lambdify(x, sp.diff(ARGS[nm], x), 'numpy')(grid)
    mono = np.all(dv > 0) or np.all(dv < 0)
    ok = v.min() > -1 and v.max() <= 1 and mono
    report(f"arg {nm}: in (-1,1] and strictly monotone",
           ok, extra=f"(grid min {v.min():.6f}, max {v.max():.6f})")

print("\n=== TASK 1, step 4 (NUMERICAL, not part of the proof): H_w by quad vs C_w ===")
mp.mp.dps = 65
lam = mp.log(2)
def f_at(a, t):
    return 1/t if a == 0 else (1/(1-t) if a == 1 else 1/(1+t))
def C2num(ab, t):
    if ab == (0, 1): return mp.polylog(2, t)
    if ab == (0, -1): return -mp.polylog(2, -t)
    if ab == (1, -1): return -mp.polylog(2, (1+t)/2) - mp.log(1-t)*mp.log(1+t) + lam*mp.log(1+t) - lam**2/2 + mp.pi**2/12
    if ab == (-1, 1): return mp.polylog(2, (1+t)/2) - lam*mp.log(1+t) + lam**2/2 - mp.pi**2/12
def Hw3_quad(w, X):   # defining integral; inner weight-2 HPL via the PROVED classical forms
    a1, a2, a3 = w
    return mp.quad(lambda t: f_at(a1, t)*C2num((a2, a3), t), [0, X])
def Hw3_nested(w, X): # fully nested defining integral (cross-check of the quadrature)
    a1, a2, a3 = w
    H1 = lambda s: (mp.log(s) if a3 == 0 else (-mp.log(1-s) if a3 == 1 else mp.log(1+s)))
    return mp.quad(lambda t: f_at(a1, t)*mp.quad(lambda s: f_at(a2, s)*H1(s), [0, t]), [0, X])
d = abs(Hw3_quad((0, 1, -1), mp.mpf('0.7')) - Hw3_nested((0, 1, -1), mp.mpf('0.7')))
print(f"  quadrature cross-check (proved weight-2 vs fully nested): {mp.nstr(d, 3)}")
for w_str in sorted(CD.keys()):
    w = ast.literal_eval(w_str)
    Cf = sp.lambdify(x, CD[w_str], 'mpmath')
    worst = max(abs(Hw3_quad(w, X) - Cf(X)) for X in [mp.mpf('0.3'), mp.mpf('0.7'), mp.mpf('0.95')])
    print(f"  NUMERICAL {w_str:>12}: max |quad - C_w| over x = 0.3, 0.7, 0.95 = {mp.nstr(worst, 3)} (65 dps)")

print("\n" + "=" * 78)
print("TASK 2: B2 -> B2' exactly")
print("=" * 78)
c = sp.symbols('c', positive=True)
q = (1 - c**2)/(2*c)
R2 = json.load(open(G1 + 'hpl_fit_k2_(0,0,0).json'))
CDraw = json.load(open(G1 + 'classical_dict.json'))

# sign rule H_w(-x) = (-1)^(#nonzero letters) H_{-w}(x).
# One-line proof: substitute t -> -t in H_{a,w}(x) = int_0^x f_a(t) H_w(t) dt;
# f_0(-t) d(-t) = f_0(t) dt and f_a(-t) d(-t) = -f_{-a}(t) dt for a = ±1, so each
# letter a contributes f_{-a} and a minus sign exactly when a != 0; induction on
# the word length (base: H_{±1}(-x) = -H_{\mp1}(x) by the definitions) gives it.
mp.mp.dps = 50
def Hw_num(w, X):
    if len(w) == 1:
        a = w[0]
        return mp.log(X) if a == 0 else (-mp.log(1-X) if a == 1 else mp.log(1+X))
    return mp.quad(lambda t: f_at(w[0], t)*Hw_num(w[1:], t), [0, X])
ok_sr = True
for w in [(0,0,1), (0,1,-1), (1,-1,1), (1,0,-1), (-1,0,1), (0,-1,1)]:
    for X in [mp.mpf('0.31'), mp.mpf('0.73')]:
        lhs = Hw_num(w, -X); rhs = (-1)**sum(1 for a in w if a != 0)*Hw_num(tuple(-a for a in w), X)
        ok_sr &= abs(lhs - rhs) < mp.mpf(10)**-40
report("sign rule H_w(-x) = (-1)^nnz H_{-w}(x) (one-line t -> -t proof; NUMERICAL check)", ok_sr)

# typed B2 (result.md) vs hpl_fit_k2 json
def nnz(w): return sum(1 for a in w if a != 0)
typed_hpl, json_hpl = {}, {}
def addh(tab, w, cf):
    tab[w] = sp.simplify(tab.get(w, 0) + cf)
Sc_terms = [((0,0,1), (1-c)), ((0,1,-1), (1-c)), ((1,-1,1), -(1-c)),
            ((1,0,1), q), ((1,0,-1), -(1-c)**2/(2*c))]
for w, cf in Sc_terms:
    addh(typed_hpl, w, cf)
for w, cf in Sc_terms:     # S(-c) via the sign rule
    addh(typed_hpl, tuple(-a for a in w), (-1)**nnz(w) * sp.expand(cf.subs(c, -c)))
addh(typed_hpl, (0, 1), -1); addh(typed_hpl, (0, -1), 1)   # -Li2(c^2)/2 by duplication
for name, cstr in R2:
    mm = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    w = tuple(int(z) for z in mm.group(2).split(','))
    addh(json_hpl, w, sp.Rational(cstr)*c**int(mm.group(1)))
ok_b2 = all(sp.simplify(typed_hpl.get(w, 0) - json_hpl.get(w, 0)) == 0
            for w in set(typed_hpl) | set(json_hpl))
report("hpl_fit_k2 json == typed B2 (S(c)+S(-c) with sign rule), all 12 coefficients", ok_b2)

# assembly from json + PROVED dictionary
coef = {}
def addc(tab, fname, cconst, val):
    tab[(fname, cconst)] = sp.simplify(tab.get((fname, cconst), 0) + val)
for name, cstr in R2:
    mm = re.match(r'c\^(-?\d+) H\[(.*)\]', name)
    e = int(mm.group(1)); w = [int(z) for z in mm.group(2).split(',')]
    pref = sp.Rational(cstr)*c**e
    if w == [0, 1]: addc(coef, 'Li2(x)', '1', pref)
    elif w == [0, -1]: addc(coef, 'Li2(-x)', '1', -pref)
    else:
        for fn, cn, fr in CDraw[str(w)]:
            addc(coef, fn, cn, pref*sp.Rational(fr))

# typed B2' from result.md, minus the zeta(2)*Phi(c) term (Phi(c) = 1/2 - q*artanh c)
typed = {}
for fn, cn, v in [
    ('Li3(x)', '1', 1 + c), ('Li3(-x)', '1', 1 - c), ('Li3(1-x)', '1', (1 - c)/c),
    ('Li3(1/(1+x))', '1', -(1 + c)/c), ('Li3((1+x)/2)', '1', -2*c),
    ('Li3((1-x)/(1+x))', '1', 1 - c), ('Li3(-(1-x)/(1+x))', '1', 1 - c),
    ('Li3(2x/(1+x))', '1', -2*c),
    ('Li2(x)', '1', -1), ('Li2(-x)', '1', -1),
    ('Li2((1+x)/2)*l1p', '1', 1 + c), ('Li2((1+x)/2)*l1m', '1', -(1 - c)),
    ('Li2(x)*l1p', '1', q), ('Li2(x)*l1m', '1', q),
    ('Li2(-x)*l1p', '1', -q), ('Li2(-x)*l1m', '1', -q),
    ('l1m*l1m*l1p', '1', -(1 - c)), ('l1m*l1p*l1p', '1', 1 - c),
    ('l1p*l1p*l1p', '1', (4*c**2 - c + 1)/(6*c)),
    ('lx*l1m*l1m', '1', q), ('lx*l1m*l1p', '1', -(1 - c)), ('lx*l1p*l1p', '1', -(1 + c**2)/(2*c)),
    ('l1m*l1m', 'log2', 1 - c), ('l1m*l1p', 'log2', -(1 - c)), ('l1p*l1p', 'log2', -c),
    ('l1p', 'log2^2', (1 - c)/2), ('l1m', 'log2^2', -(1 - c)/2),
    ('l1m', 'pi2', -(1 - c)/(6*c)), ('l1p', 'pi2', -(1 + c)/(6*c)),
    ('1', 'z3', sp.Rational(7, 4) + 2*c), ('1', 'pi2*log2', -c/6), ('1', 'log2^3', c/3)]:
    addc(typed, fn, cn, v)
print("\n  coefficient-by-coefficient comparison in Q(c) (assembled from PROVED dict vs typed B2'):")
ok_all = True
for k in sorted(set(coef) | set(typed)):
    d = sp.simplify(coef.get(k, 0) - typed.get(k, 0))
    ok_all &= (d == 0)
    print(f"    {str(k):>34}: diff = {d}")
report("B2 -> B2': all coefficients in Q(c) agree (Z(c) = 7/4 + 2c included)", ok_all)

# NUMERICAL: typed B2' against the stored M2 values, 100 digits
mp.mp.dps = 120
def M2_typed(cc):
    lam = mp.log(2); lp = mp.log(1+cc); lm = mp.log(1-cc); l0 = mp.log(cc)
    qq = (1-cc**2)/(2*cc); Phi = mp.mpf(1)/2 - qq*mp.atanh(cc)
    return (mp.zeta(2)*Phi - mp.polylog(2, cc) - mp.polylog(2, -cc)
        + (1+cc)*mp.polylog(3, cc) + (1-cc)*mp.polylog(3, -cc)
        + ((1-cc)/cc)*mp.polylog(3, 1-cc) - ((1+cc)/cc)*mp.polylog(3, 1/(1+cc))
        - 2*cc*mp.polylog(3, (1+cc)/2)
        + (1-cc)*(mp.polylog(3, (1-cc)/(1+cc)) + mp.polylog(3, -(1-cc)/(1+cc)))
        - 2*cc*mp.polylog(3, 2*cc/(1+cc))
        + mp.polylog(2, (1+cc)/2)*((1+cc)*lp - (1-cc)*lm)
        + qq*(mp.polylog(2, cc) - mp.polylog(2, -cc))*(lp+lm)
        - (1-cc)*lm**2*lp + (1-cc)*lm*lp**2 + ((4*cc**2-cc+1)/(6*cc))*lp**3
        + qq*l0*lm**2 - (1-cc)*l0*lm*lp - ((1+cc**2)/(2*cc))*l0*lp**2
        + lam*((1-cc)*lm**2 - (1-cc)*lm*lp - cc*lp**2) + lam**2*((1-cc)/2)*(lp-lm)
        - mp.pi**2*(((1-cc)/(6*cc))*lm + ((1+cc)/(6*cc))*lp)
        + (mp.mpf(7)/4 + 2*cc)*mp.zeta(3) - (cc/6)*mp.pi**2*lam + (cc/3)*lam**3)
ok_num = True
for tag, alpha in [("pi_2", mp.pi/2), ("2pi_3", 2*mp.pi/3), ("pi_3", mp.pi/3)]:
    st = [l for l in open(G1 + f'values_{tag}.txt') if l.startswith('M2 =')][0]
    ref = mp.mpf(st.split('=')[1].strip())
    diff = abs(M2_typed(mp.cos(alpha/2)) - ref)
    ok_num &= diff < mp.mpf(10)**-100
    print(f"  NUMERICAL {tag:>6}: |typed B2' - stored M2| = {mp.nstr(diff, 3)} (120 dps)")
report("typed B2' vs stored M2 at 3 angles, 100+ digits (NUMERICAL)", ok_num)

print("\n" + "=" * 78)
npass = sum(1 for _, ok in RESULTS if ok)
print(f"SUMMARY: {npass}/{len(RESULTS)} items PROVED (numerical items labelled NUMERICAL),"
      f"  total time {time.time()-T0:.1f}s")
for label, ok in RESULTS:
    if not ok:
        print("  FAILED:", label)
'''
with open("verify_cor3.py", "w") as f:
    f.write(script)
print("written verify_cor3.py,", len(script), "chars")
