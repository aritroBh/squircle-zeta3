import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
summary = """# AV5b final summary

## Layer Lemma (PROVED)
h_tau = m·A·exp(Σ tau^k g_k(s)), s = (theta-theta0)/tau, tent m = H cos u + T|sin u| (alpha = T/H),
bulk A = exp(tau a1 + tau^2 a2 + ...), g_k decaying. P = int h dtheta:
c2_lay = H int g1 ; c3_lay = H int [g2 + g1^2/2 + (alpha|s| + a1) g1] ; c4_lay = H int [G3 + alpha|s|G2 - s^2 g1/2].
tau = 1-1/p = d-d^2+d^3-... ; P'''(1) = 6(c1 - 2 c2 + c3).

## Softplus moments (PROVED)
int_0^inf s^j g1 = j! eta(j+2)/2^{j+1}; int_0^inf g1^2 = zeta(3)/8;
int_0^inf s^a g1^b = a!/2^{a+1} * alternating Euler sum (coloured MZV, weight a+b+1, depth b).
J = int_0^inf (s g1^2/2 + g1^3/6) = zeta(4)/32 (NUMERICAL, 60 digits, resid -4.6e-60).

## Criterion (proved if / proved for subclasses only-if; see result.md)
softplus g1 (moments are eta values) + non-cancelling combination <=> zeta(3) at tau^3.

## Paper checks (own code)
- m=2: c2 = sqrt2 z2, c3 = sqrt2 z3 (PROVED); P'''(1) = 6 sqrt2 z3 - 2 sqrt2 pi^2 (jet, 55 digits).
- m=3: g1 = ln2 + (1/2) log(1+e^{-4 alpha |s|}) (derived, 1e-12); int g2 = -sqrt3 z3/8 (60 digits);
  c3_lay-zeta = -sqrt3 z3/32 = C3 z3/36; C3 = -9 sqrt3/8 reproduced. P'(1) = 2 sqrt3 ln2 exact;
  P''(1), P'''(1) by PSLQ at 95 digits with C3 as zeta(3) coefficient.
- m=4,5,6: P'(1) = paper R_1 (42-45 digits); P'''(1): m=4 -4.39449630527605215997911011380609991382237745,
  m=5 -4.77620018088223841663202245344488583697077, m=6 -5.20324042951766621029582970427243461666281
  (= paper to 21 digits). C_m formula verified as values.

## Task C
- F1: h = m*(sum_v e^{(v.n-m)/tau})^tau (square): P = 8 + (pi^2/3) tau^2 + 2 zeta(3) tau^3 + (zeta(4)/4) tau^4
  (c3 to 21 digits, c4 to 14 digits; convexity min h+h'' ~ 0 at dps 55).
- F2: h = tau log sum_v e^{v.n/tau}: P = 8 + (pi^2/3) tau^2 + 0 tau^3 + (7 zeta(4)/8) tau^4 (c3 = 0 to 9e-11).

## verify_av5b.py: 18/18 PASS in 22s.
"""
with open("summary.md","w") as f: f.write(summary)
print("summary.md updated")
print(open("progress.log").read())
