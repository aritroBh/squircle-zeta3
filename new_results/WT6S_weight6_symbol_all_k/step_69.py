
import os
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
md = r'''# WT6S: the weight-six part W6(k) of the squircle coefficients M_k — E_{5,1..5} for 6 ≤ k ≤ 52

**Verdict:** W6(k) = Σ_{m=1}^{5} C(k−6, m−1)·E_{5,m} holds exactly for 6 ≤ k ≤ 52, with E_{5,m} fixed
Q(√2)-combinations of the twelve integrals N_{(i,j),r} = ∫_b^a log^i t · log^j(1−t) dv/(v−r)
(i+j=5, r ∈ {0, a=1+√2, −a, −b=1−√2}). **E_{5,5} and E_{5,4} are classical at the symbol level
(all (3,3)/(4,2) pairings 0); E_{5,1}, E_{5,2}, E_{5,3} are NOT classical**, with explicit depth-two
Li_{3,3}/Li_{4,2} corrections on M5b's four-atom argument set. The depth-two part of W6(k) is thus a
fixed four-atom combination with coefficients polynomial in k of degree ≤ 2 (one atom, Li_{4,2}(b,1/√2),
has the constant coefficient 120√2 for all k). E_{5,5} = 4!·c(6), confirmed to 151 digits.
All controls pass; see the evidence table. Nothing below is claimed that my own code did not print.

Notation: λ = log 2, L = log(1+√2), X = λ-letter, Y = L-letter, P = [X,Y], R = √2, a = 1+√2, b = √2−1,
t = (v²+2v−1)/(4v). Pairings: (3,3) against [X,P]∧[Y,P]; (4,2) = (a,b,c) against {[X,[X,P]],[Y,[X,P]],[Y,[Y,P]]}∧P
(same normalisation as M5a). "Classical at the symbol level" = vanishing (3,3) and (4,2) Lie pairings
(necessary, not sufficient; it is the depth-two test).

## A. E_{5,m} from the kernels (EXACT, k = 6..52)

The degree-5 kernel residues res_{ij,r}(k) = rat(k) + √2·sq(k) (quartic in k, TOWN1
`deg5_residue_polynomials.json`) were rewritten in the binomial basis C(k−6, m−1) exactly
(`step_03.py`): res_{ij,r}(k) = Σ_m C(k−6,m−1)·e_{ij,r,m}, so E_{5,m} = Σ_{ij,r} e_{ij,r,m}·N_{(i,j),r}.
The rewrite reproduces the residue polynomials **exactly** at k = 6, 7, 10, 21, 52 (A1, PASS), and
e_{ij,r,5} = 24·(k⁴-coefficient of res_{ij,r}) for every cell (PASS). Full table:
`e5m_binomial_coeffs.json`. The two ends:

- **E_{5,1}** = (1/640+5R/384)N(5,0,0) + (5/128)N(5,0,a) − (193/1920)N(5,0,−a) + (7/120)N(5,0,−b)
  − (13/384+23R/384)N(4,1,0) − (97/384)N(4,1,a) + (187/384)N(4,1,−a) − (1/6)N(4,1,−b)
  + (161/384+3R/64)N(3,2,0) + (107/384)N(3,2,a) − (161/192)N(3,2,−a) − (107/384)N(3,2,−b)
- **E_{5,5}** = (2/15)N(5,0,a) − (2/15)N(5,0,−a) − (2/3)N(4,1,a) + (2/3)N(4,1,−a)
  + (2/3)N(3,2,0) + (2/3)N(3,2,a) − (4/3)N(3,2,−a) − (2/3)N(3,2,−b)   (8 terms, all rational)
- (E_{5,2}, E_{5,3}, E_{5,4}: 12 terms each; see `e5m_binomial_coeffs.json`.)

## B. Symbols and cobracket pairings (EXACT; controls first)

The weight-6 symbol of each N_{(i,j),r} was computed exactly (`step_04.py`, `n6_symbols_pairings.json`;
29–64 words each over the alphabet {X, Y}). Pairings are Q(√2)-linear, so pv(E_{5,m}) = Σ e_{ij,r,m}·pv(N).

**CONTROLS (both PASS):**
- **(i)** W6 at k=5 built from `kern_5.pkl` (degree-5 keys) gives (3,3) = **449√2/4**, (4,2) = **(−30, 15√2, −9)**
  — exactly M5a's values.
- **(ii)** At weight five (deg4 kernels in the binomial basis at k=5; the e's reproduce WT5's
  `deltas_deg4.json` exactly), the (3,2) pairings are E1 = **(−17/2, 0)**, E2 = **(−2, 0)**, E3 = **(0, 0)**
  — exactly WT5's.

**Pairings of the E_{5,m} (33, 42a, 42b, 42c):**

| E_{5,m} | (3,3) | (4,2)a | (4,2)b | (4,2)c | classical? |
|---|---|---|---|---|---|
| E_{5,1} | 46√2 | −75/2 | 15√2 | 0 | NO |
| E_{5,2} | −21√2/2 | −35/2 | 0 | 21 | NO |
| E_{5,3} | −6√2 | −10 | 0 | 12 | NO |
| E_{5,4} | 0 | 0 | 0 | 0 | YES (symbol level) |
| E_{5,5} | 0 | 0 | 0 | 0 | YES (symbol level) |

## C. Which E_{5,m} are classical; pairings of W6(k) as polynomials in k (EXACT)

**E_{5,5} and E_{5,4} are classical at the symbol level** (all four pairings exactly 0);
E_{5,1}, E_{5,2}, E_{5,3} are not. (The task expected E_{5,5}; E_{5,4} being classical too is stronger,
and it lowers the depth-two degree.) Since E_{5,4} = E_{5,5} = 0 in the pairings, the W6(k) pairing
vector is a polynomial in k of degree **≤ 2** (not merely ≤ 3). Exact power-basis forms (`step_05.py`),
each verified against the full weight-6 symbol built independently from the kernels at k = 6, 7, 10, 21, 52:

- P_{3,3}(k) = (−3k² + 57k/2 − 17)·√2
- P_{4,2a}(k) = −5k² + 95k/2 − 285/2
- P_{4,2b}(k) = 15·√2   **(constant in k)**
- P_{4,2c}(k) = 6k² − 57k + 126

*Curiosity (EXACT):* extrapolating the quartic family to k = 5 gives (−30, 15√2, −9) in the three (4,2)
components — matching the actual kern_5.pkl/N6 values — but (3,3) = 101√2/2 ≠ 449√2/4. So the true k=5
kernels leave the quartic family only in the (3,3) direction (difference −247√2/4). (`step_09.py`.)

## D. Depth-two parts N(E_{5,m}) (EXACT; unique)

Atom pairing vectors computed independently (`step_06.py`) match M5a's table, and M5b's N6 combination
reproduces the k=5 target (control, PASS). The pairing matrix of the four atoms has det = **9/512 ≠ 0**,
so **every** target has a **unique** depth-two representative on this set (no enlargement needed):

- **N(E_{5,1})** = (−56√2)·Li₃,₃(b,1) + (−200+80√2)·Li₃,₃(1/√2,1) + (16√2)·Li₃,₃(b,2) + (120√2)·Li₄,₂(b,1/√2)
- **N(E_{5,2})** = (28−42√2)·Li₃,₃(b,1) − (280/3)·Li₃,₃(1/√2,1) − (21√2/2)·Li₃,₃(b,2)
- **N(E_{5,3})** = (16−24√2)·Li₃,₃(b,1) − (160/3)·Li₃,₃(1/√2,1) − (6√2)·Li₃,₃(b,2)
- **N(E_{5,4}) = 0, N(E_{5,5}) = 0**

Each verified by direct recomputation of the combination's pairings (PASS). Depth-two of W6(k)
(`step_09.py`, `n_w6k_polynomial.json`; pairings match P(k) at k = 6, 7, 10, 21, 52):

- Li₃,₃(b,1):    (8k²−76k+168) + (−12k²+114k−308)·√2
- Li₃,₃(1/√2,1): (−80k²/3+760k/3−760) + 80·√2
- Li₃,₃(b,2):    (−3k²+57k/2−47)·√2
- Li₄,₂(b,1/√2): **120·√2  (constant in k)**  ← sole source of the constant P_{4,2b} = 15√2

## E. Numerics (150 digits; NUMERICAL)

The twelve N_{(i,j),r} were integrated by two independent methods (`step_07.py`): raw tanh-sinh over
[b, 1/2, 1, 2, a] (M5a setup) and a power-substitution method (split [b,1]+[1,a], v = b+(1−b)u¹² /
v = a−(a−1)u¹², removing the endpoint log singularities; for r = a on [1,a] the pole is removed
analytically via dv/(v−a) = −(12/u)du). Worst |raw − sub| = 1.65e−162 over the 12 integrals; stored in
`n6_values.json`. E_{5,m} = Σ e·N assembled (`step_08.py`, `e5m_values.json`, 150 digits):

- E_{5,1} = −4.418706461977750371981415046115271660794369502753030878…
- E_{5,2} = −6.292447336030847786607049789803405607923022406569969455…
- E_{5,3} =  7.338991964028723995633475132236234412143401235900679477…
- E_{5,4} = 17.957475903672145702821029549088656826649556951836575265…
- E_{5,5} =  8.056501397376488568184134716244283853076847305693403011…

Checks: **c(6) = E_{5,5}/24 matches LEAD1's 405-digit value to 2.3e−151**; W6(6), W6(7), W6(10) by direct
quadrature of the full kern_k.pkl degree-5 kernels equal Σ_m C(k−6,m−1)E_{5,m} to ~1e−148.
The classical remainders of E_{5,1..4} were NOT attempted (per the task); E_{5,5} = 24·c(6) is the
LEAD-law Γ-world constant.

## F. Deliverables and rerun

`verify_wt6s.py` reruns A–E in one process (shared symbol cache) in ~5 min, printing PASS/FAIL per item
(last run: ALL PASS, 293 s). Steps: `step_03.py` (A), `step_04.py` (B + controls), `step_05.py` (C),
`step_06.py` (D), `step_07.py`/`step_08.py` (E), `step_09.py` (D-poly + k=5 observation). Core module
`wt6s_core.py` (adapted from M5a `w6_symbol.py`, same normalisation); copied (never edited elsewhere):
`qsqrt2.py`, `w4_symbol.py`, `deg5/deg4_residue_polynomials.json`, `kern_{5,6,7,10}.pkl`,
`kernels_deg4.json`, `deltas_deg4.json`, `omega5.json`, `c6_lead_reference.json` (LEAD1 405-digit c(6)).
Data: `e5m_binomial_coeffs.json`, `e5m_lead_coeffs.json`, `n6_symbols_pairings.json`, `e5m_pairings.json`,
`w6k_pairing_polynomials.json`, `n_e5m_solution.json`, `n_w6k_polynomial.json`, `n6_values.json`,
`e5m_values.json`. Rerun:

```
PY='<project>/.venv/bin/python'
nice -n 10 $PY -B verify_wt6s.py        # all controls + A-E, ~5 min
```

## Evidence table

| item | statement | type | result |
|---|---|---|---|
| A1 | binomial rewrite of 15 residue cells == polynomials, k=6,7,10,21,52 | EXACT | PASS |
| A2 | e_{ij,r,5} = 24·k⁴-coeff of res_{ij,r}, all 15 cells | EXACT | PASS |
| B-(i) | W6(k=5) from kern_5.pkl pairings | EXACT | (449√2/4, −30, 15√2, −9) = M5a |
| B-(ii) | weight-5 E1,E2,E3 (3,2) pairings | EXACT | (−17/2,0),(−2,0),(0,0) = WT5 |
| B | 12 N_{(i,j),r} pairing vectors + E_{5,m} pairings | EXACT | table above |
| C | E_{5,4}, E_{5,5} all pairings 0 (classical, symbol level) | EXACT | PASS |
| C | P(k) == independent full symbol at k=6,7,10,21,52 | EXACT | PASS (all 5 k) |
| D | atom pairing vectors == M5a table; M5b N6 control | EXACT | PASS |
| D | det(pairing matrix) = 9/512 ≠ 0 (unique N) | EXACT | PASS |
| D | N(E_{5,m}) reproduces all pairings, m=1..5 | EXACT | PASS ×5 |
| E | N values raw vs sub, worst over 12 | NUMERICAL | 1.65e−162 |
| E | c(6) = E_{5,5}/24 vs LEAD1 405-digit | NUMERICAL | 2.3e−151 (151 digits) |
| E | W6(6/7/10) direct kern_k.pkl vs Σ C(k−6,m−1)E_{5,m} | NUMERICAL | ~1e−148 |

**Bottom line:** For 6 ≤ k ≤ 52, W6(k) is a quartic in k whose two top coefficients E_{5,5} (= 24c(6),
Γ-world) and E_{5,4} are classical at the symbol level, while E_{5,1}, E_{5,2}, E_{5,3} carry explicit,
unique depth-two Li₃,₃/Li₄,₂ corrections on M5b's four atoms with polynomial (degree ≤ 2) coefficients —
Li₄,₂(b,1/√2) entering with the constant coefficient 120√2. The k = 5 point (M5a/M5b's N6) sits off the
quartic family only in its (3,3) pairing.
'''
open(os.path.join(WD, 'result.md'), 'w').write(md)
print(f"result.md written, {len(md.splitlines())} lines")
