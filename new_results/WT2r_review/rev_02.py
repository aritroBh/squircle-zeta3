
import sys, pickle
sys.path.insert(0,'.')
from fastred import *
import sympy as sp
from fastred import v
I = sp.I; S2 = sp.sqrt(2)

# ---- Lemma A ingredient 1: DL1, DL2 zeros at +-i and derivatives (EXACT, symbolic) ----
a = 1+S2; b = S2-1
DL1s = 1/(v-b) + 1/(v+a) - 1/v
DL2s = 1/(v-a) + 1/(v+b) - 1/v
for name, DL in [('DL1',DL1s), ('DL2',DL2s)]:
    at_i  = sp.simplify(DL.subs(v, I)); at_mi = sp.simplify(DL.subs(v, -I))
    d_i   = sp.simplify(sp.diff(DL,v).subs(v, I)); d_mi = sp.simplify(sp.diff(DL,v).subs(v, -I))
    print(f"{name}: value at +i = {at_i}, at -i = {at_mi};  deriv at +i = {d_i}, at -i = {d_mi}")
print("tau1_0 claim (-1-i)/2:", sp.simplify(sp.diff(DL1s,v).subs(v,I) - (-1-I)/2)==0)
print("tau2_0 claim (-1+i)/2:", sp.simplify(sp.diff(DL2s,v).subs(v,I) - (-1+I)/2)==0)
print("DL1'(-i) == conj(DL1'(i)) == tau2_0:", sp.simplify(sp.diff(DL1s,v).subs(v,-I) - (-1+I)/2)==0)
print("DL2'(-i) == tau1_0:", sp.simplify(sp.diff(DL2s,v).subs(v,-I) - (-1-I)/2)==0)

# ---- Lemma A ingredient 2: F preserves top pole order with factor -1 (exact, m=1..40) ----
# F: R(v) -> R(1/v)/v^2.  R = c/(v+i)^m  =>  F(R) = -c/(v-i)^m + lower order terms at v=i.
u = sp.Symbol('u')
okF = True
for m in range(1,41):
    c = sp.Symbol('c')
    R = c/(v+I)**m
    FR = sp.together(R.subs(v, 1/v)/v**2)
    # series at v = I: substitute v = I + u, expand in u
    ser = sp.series(FR.subs(v, I+u), u, 0, -m+2).removeO().expand()
    # coefficient of u^(-m)
    coeff_top = sp.expand(ser).coeff(u, -m)
    # remainder should have only powers u^(-m+1) and higher
    rem = sp.expand(ser - coeff_top*u**(-m))
    powers = [t.as_base_exp()[1] for t in sp.Add.make_args(rem)] if rem != 0 else []
    minpow = min([sp.degree(p, u) if p.is_Pow else (1 if p==u else 0) for p in powers], default=0)
    okF &= (sp.simplify(coeff_top + c)==0) and all((p.is_Pow and p.exp>=-m+1 and p.base==u) or p==u or (not p.has(u)) for p in powers)
    if m<=3: print(f"  m={m}: F(c/(v+i)^{m}) top coeff at +i = {coeff_top} (expect -c), rest min power ok: {okF}")
print("F: c/(v+i)^m -> -c/(v-i)^m + O((v-i)^(-m+1)) for m=1..40:", okF)

# ---- Lemma A ingredient 3: A/B exponential blocks have NO +-i poles, for all k (structural) ----
# AB_part(v) = -(v^2+1)/(2 k! v (v-b)(v+a)): denominator v(v-b)(v+a) -- real poles only. Symbolic for k=1..6:
import buildintegrand
okAB = True
for k in range(1,7):
    P = buildintegrand.AB_part_PF(k)
    okAB &= (not P[1]) and (not P[2])   # empty ipos/ineg lists
print("AB_part_PF(k) has empty +-i principal parts for k=1..6:", okAB, " (denominator v(v-b)(v+a): real poles only, all k)")

# ---- Lemma A ingredient 4: sqrt block base-case order <= 2d-2 (structural from build) ----
# sqrt_part_PF denominator = (v^2-a^2)(v^2-b^2) v^e (v^2+1)^{2d-2}: so +-i order <= 2d-2. Verify PF structure d=2..9:
from genintegrand import S_poly
okSQ = True
for d in range(2,10):
    for m in range(0,d+1):
        P = buildintegrand.sqrt_part_PF(d, S_poly(d,m))
        if len(P[1])>2*d-2 or len(P[2])>2*d-2: okSQ=False
print("sqrt_part_PF(d,m): +-i order <= 2d-2 for d=2..9 all m:", okSQ)
