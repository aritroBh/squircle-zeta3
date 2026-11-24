res = """# D13t — adversarial review of the D13s fixes (Kimi K3) — WORKING DRAFT

Theorem under review: for all N >= 100, |M_k^main + (sqrt(pi)/N) e^{-sqrt(pi N)} cos sqrt(pi N)| <= C N^{-3/2} e^{-sqrt(pi N)}, C = 254.2276 (certified sup 254.2275176).

## Verdicts (updated as items complete)

| item | content | verdict |
|---|---|---|
| R1 | B2 fix: A1 enclosure in fixed/step_38.py | **HOLDS** (proof + 110884-pt containment) |
| R2 | B1 fix: other rt2-for-2 slips | **HOLDS** (grep+diff+sympy audit; no other slip) |
| R3 | H1: B_true exact regular form | pending |
| R4 | H2: |sH_c(s)|<=1 proof | pending |
| R5 | H3: certify_H3 box coverage, t<=0.12731 | pending |
| R6 | rerun certificates | pending |

## R1 (B2 fix) — HOLDS

- sympy independently: f(x)=1-atan(x)/x, f'=g/x^2, g=atan(x)-x/(1+x^2), g(0)=0, g'=2x^2/(1+x^2)^2 >= 0 (EXACT, residual 0) => f nondecreasing on [0,oo); the enclosure A1 = f(zz)/y in [f(zz_lo)/y_up, f(zz_up)/y_lo] is a valid outer enclosure (arb boxes outer-enclose the range; dependency only loosens).
- 110884 sample (y,sigma) points over 149 cells x 4 sigma-subcells (N near 100, 1891 [sweep-worst cell 1890.7128907662395,1891.0910333443928], 1e4, 1e6, 1e8, 1e9, 1e10-eps; wide and ratio-1.0002 boxes): worst A1 containment violation = -2.09e-9 (first batch) / -2.30e-10 (random cells) => CONTAINED everywhere. zz-box worst violation -1.07e-12.
- saddle sigma* in [0.42858, 0.5] for N in [100,1e20] => the true saddle always lies inside the sigma-box [0,1].
- phi = (pi/rt2) B_h[A1'/y] independently confirmed (my own Phi, saddle root-finding): worst residual 4.44e-56 @60 dps over 2006 points incl N=100 and N=1891. This is exactly the phase coded (T2n = 2y^2 sg - rt2 rho y - 2 rho A1m y; Bh = -y sg/(2rho) + T2n/(2y(2+rt2 rho y)); phireg = (pi/rt2) Bh).
- Note: a first version of MY test produced a spurious +0.609 violation — caused by my own inverted y-box (r<1), not by the certificate. Redone with valid boxes: all contained.

## R2 (B1 fix) — HOLDS

- grep `rt2/`, `.sqrt()`, `**arb`, `ecw`, `c_w` over fixed/step_10_J192.py and fixed/step_38.py; diff vs buggy baselines: ONLY the intended fix lines changed (B1 ecw; B2 A1m; H1 breg->Btrue box; H2 TL0; TLfN (pi/2)^1.5 factor). step_18.py identical.
- ecw fixed = 0.483941449038287 = sqrt(2/pi) e^{-1/2} CORRECT (buggy 0.406944629688983; ratio = 2^{1/4} = 1.189207115 exactly, matching D13s's table).
- Every other audited constant: GRreg 2^{9/4} = (2rt2)^{3/2} OK; GLm2reg 2^{19/4} = 2^{9/4}·2^{5/2} OK (eta>=1/2 factor); G0reg = sqpi/4·2^{3/2} OK; mSreg factor (N y0^3/2)/y = pi rho^3/4 OK; TL0 exponent form identical to 2(2+rt2)^{-N/2}e^{sqrt(pi N)}/(N-2); d6 = (pi/rt2)·0.16666814 = 0.3702435178 OK.
- sympy EXACT: (t^2+rt2 t)^2+(rt2 t+1)^2 = (t^2+rt2 t+1)^2 (so |acb(t^2+rt2 t, rt2 t+1)| = t^2+rt2 t+1, p2mag line correct); H2 numerator 3x^2+6rt2 kx+4k^2 > 0; AM-GM residual (x-1)^2 >= 0.
"""
open('result.md','w').write(res)
print("result.md written,", len(res), "chars")
