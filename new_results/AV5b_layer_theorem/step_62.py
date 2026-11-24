import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
summary = """# AV5b: Layer theorem — running summary

## Task A (Layer Lemma + softplus moments)
- Setting: h_tau(theta) = m(theta) * A(theta,tau) * exp(sum tau^k g_k((theta-theta0)/tau)), m = polygon tent
  (H = height, slopes +/- T, alpha = T/H), A = bulk analytic factor exp(tau a1 + tau^2 a2 + ...), g_k decaying layer profiles.
- Perimeter P = int h dtheta (Cauchy). Expansion: P = P_poly + c1 tau + c2 tau^2 + c3 tau^3 + ...
  c2_layer = H int g1 ds ; c3_layer = H int [g2 + g1^2/2 + (alpha|s| + a1) g1] ds  [VALIDATED m=2 to 7 digits, m=3 below]
  c4_layer = H int [G3 + alpha|s|G2 - s^2 g1/2] ds, G3 = g3 + g1 g2 + g1^3/6.
- Softplus g1 = log(1+e^{-2|s|}): int_0^inf s^j g1 ds = j! eta(j+2)/2^{j+1} (PROVED: geometric series + Fubini);
  int_0^inf g1^2 ds = zeta(3)/8 (PROVED: u=e^{-2s}, int_0^1 ln^2(1+u)/u du = zeta(3)/4 via Euler sum 5z3/8);
  int_0^inf s^a g1^b ds = a!/2^{a+1} * sum_{k1..kb>=1} prod(-1)^{ki+1}/ki / (k1+...+kb)^{a+1}
  = alternating (coloured, 2nd roots of unity) MZV of weight a+b+1, depth b (PROVED, abs. conv. by AM-GM).

## Task B (checks against the paper)
- m=2 (squircle/diamond): c2 = sqrt2 zeta(2), c3 = sqrt2 zeta(3) (60 digits, matches AV5/paper pure-weight);
  P'''(1) = 6 sqrt2 zeta(3) - 2 sqrt2 pi^2 — my jet engine reproduces to 55 digits [step_14/18].
  tau = 1-1/p = d-d^2+d^3-..., d = p-1: P'''(1) = 6(c1 - 2 c2 + c3) = 6A3, A3 = coeff of d^3.
- m=3 (hexagon, fully local => local = total):
  P(1) = 2 sqrt3 (exact), P'(1) = 2 sqrt3 ln 2 (jet engine, matches paper R_1 to 60 digits);
  P''(1) = sqrt3(pi^2 + 16 ln^2 2 - 32 ln 2)/8 [PSLQ, 98 digits];
  P'''(1) = -9 sqrt3/8 zeta(3) + sqrt3(3 pi^2 ln2/8 + 2 ln^3 2 - 3 pi^2/4 - 12 ln^2 2 + 12 ln 2) [PSLQ, 95 digits]
  => zeta(3) coefficient = C_3 = -9 sqrt3/8 EXACTLY (paper's Thm C_m, fully local).
  Layer structure: bulk a1 = ln2 (const on fan), a2 = 0;
  g1(s) = ln2 + (1/2) log(1+e^{-4 alpha |s|}), alpha = tan(pi/6) = 1/sqrt3 — DERIVED analytically from the
  two-kink layer equation t* = sigma(4 alpha s) (PROVED) and verified pointwise to 1e-12 [ckpt_prof_m3].
  c3_lay per corner = (1/2)[int g2 + sqrt3 z3/16 + sqrt3 pi^2 ln2/48]; backed-out int g2 = -sqrt3 zeta(3)/8
  to ALL 60 digits (from exact A3 = P'''/6 and closed-form pieces) [av5b_g2back].
  => c3_lay-zeta = (1/2)[-sqrt3/8 + sqrt3/64 + 3 sqrt3/64] zeta(3) = -sqrt3 zeta(3)/32 = C_3 zeta(3)/36.
  => the Layer Lemma's corner computation reproduces the paper's C_3 EXACTLY. Sign: negative comes from
  the g2 profile moment (-sqrt3 z3/8) beating the g1^2 and |s|g1 moments (+4 sqrt3 z3/64).
- m=4,5,6: paper totals contain non-local Li_3 parts; local part = C_m (basis-independent).
  m=4: paper total zeta(3) coeff = (87/2 - 117 sqrt2/4) u, u = sqrt(2+sqrt2); local C_4 = (54-39 sqrt2) u.
  m=6: P'''(1) = -5.20324042951766621 (paper, 630-digit quadrature).
  Profiles extracted [ckpt_prof_m4]: a1 = log(1+sqrt2) + sqrt2 ln2/(2(1+sqrt2)) (vertex value, matches 16 digits),
  a2 = -0.014572182, a3 = -0.014860409, D1(0) = 1.297869118.

## Engines (all in workdir)
- Jet engine (exact p-derivatives of polar speed integrand at p=1 via order-3 jets): P^{(k)}(1), validated
  m=2 (55 digits vs paper), m=3 (P' exact match). [step_18]
- Richardson extraction of DeltaP(d) series. [step_30, av5b_fit_true]
- Profile tower extraction D1,D2,D3(s) from support function. [av5b_profiles.py, ckpt_prof_m*.json]
- Arc-length two-kink local model (exact for m=2, 35 digits vs squircle). [step_24/25]

## TODO
- Task C: two new families (multiplicative softplus square -> c3 = 2 zeta(3); vertex-LSE square -> c3 = 0).
- m=4..6 totals via jet engine; m=4 corner back-out int g2 -> C_4/48 check.
- verify_av5b.py, result.md, FINAL report.
"""
with open("summary.md","w") as f: f.write(summary)

# Jet engine P'''(1) for m=4
src = open("step_18.py").read()
body = src.split("print(")[0]
step = body + '''
for k in [1,2,3]:
    v = P_deriv(4, k)
    print("m=4 P^(%d)(1) = %s" % (k, mp.nstr(v, 45)))
'''
exec(step)
