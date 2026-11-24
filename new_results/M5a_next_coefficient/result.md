# M5a — the weight-6 part of M5 needs depth two (symbol level), exactly as W5 did

**Verdict: W6 is NOT classical.** Its exact symbol obstruction at weight 6 is nonzero in both the (3,3) and (4,2) Lie parts:
(3,3) pairing against [X,P]∧[Y,P] = **449√2/4**; (4,2) pairings against B4⊗P = **(−30, 15√2, −9)** for B4 = {[X,[X,P]], [Y,[X,P]], [Y,[Y,P]]} (X = λ = log 2, Y = L = log(1+√2), P = [X,Y]). All controls pass. A depth-two correction
**N6 = (−12+209√2)·Li₃,₃(b,1) + (−160+80√2)·Li₃,₃(1/√2,1) + (329√2/4)·Li₃,₃(b,2) + 120√2·Li₄,₂(b,1/√2)**
reproduces all four pairings exactly (symbol-level analogue of M4c's N5). So M6's top-weight part needs Li_{4,2}/Li_{3,3}, as M4's needed Li_{3,2}.

Everything below was printed by my own code (work dir `new_results/M5a_next_coefficient/`).

## 1. Exact Ω's for k = 5 (`ibp_reduce.py 5`, 27 s; parsed to `omega5.json`, exact reconstruction asserted)

M5 = Σ ∫_b^a l1ⁱ l2ʲ Ω_ij(v) dv, a = 1+√2, b = √2−1:

- **W6** (i+j=5): Ω₃₂ = (7√2/384)/v; Ω₄₁ = (−17√2/768)/v; Ω₅₀ = (1/120 − √2/3840)/v + (−1/60)/(v+a).
- **W5″** (i+j=4): Ω₂₂ = (27/64 + 9√2/64)/v + (−27/32)/(v+a); Ω₃₁ = (7/48 − 7√2/48)/v + (5/48)/(v−a) + (13/48)/(v+a) + (−2/3)/(v+b); Ω₄₀ = (−5/192 + √2/192)/v + (−11/192)/(v−a) + (−1/64)/(v+a) + (1/8)/(v+b).
- **W4‴** (i+j=3): Ω₂₁ = (23/16 + √2/8)/v + (17/16)/(v−a) + (−23/8)/(v+a) + (−17/16)/(v+b); Ω₃₀ = (−1/24 − √2/8)/v + (−5/8)/(v−a) + (25/24)/(v+a) + (−1/3)/(v+b).
- **W3⁽⁵⁾** (i+j=2): Ω₁₁ = (5/4 + √2/4)/v + (−5/2)/(v+a); Ω₂₀ = (1/4 − √2/4)/v + (1/4)/(v−a) + (3/4)/(v+a) + (−3/2)/(v+b).
- **T** (i+j=1): Ω₁₀ = (1/2)/v + (1/2)/(v−a) + (−1)/(v+a) + (−1/2)/(v+b) — **identical to the k=4 Ω₁₀**.

## 2. Weight split and checks (quadrature at 100 dps, err est ~1e-146)

M5 = W6 + W5″ + W4‴ + W3⁽⁵⁾ + T:

| piece | value (90 digits) |
|---|---|
| W6 | −0.907608581655986804238627517183657700901005366247507315345945374787229055230616736626224053 |
| W5″ | 1.73554296628742416074075754504308245362680589046550042341518215466981343227993905904899127 |
| W4‴ | 0.874350842195575664885373894182087417876473168463622868304740230982501749385502482900429188 |
| W3⁽⁵⁾ | −2.2304703483082077291888754667732304301926353433343586141458505079006376155076385467009693 |
| T | 0.531346770191606965284879423243172540237375713200336006150605590748418213165758836508394291 |
| M5 | 0.00316164871041225748350787851145428064701406254759336837873209371286672409294509513062140084 |

- pieces(quadrature) − M5(full k=5 integrand) = **1.0e-126**; with exact tails instead = **4.0e-100** (110-dps check).
- **Neville extrapolation** of the perimeter (independent method; using proved M0,M1, paper's 1200-digit M2,M3, M4a's 800-digit M4): M5 Neville − M5 quad = **2.2e-56**.
- **Weight-2 piece: YES, T = (π² + 6λ²)/24 again** — exact symbolic zero in the tails representation (`T - (pi^2+6log^2 2)/24 = 0`) and 0.0 at 110 dps.
- **W3⁽⁵⁾: YES, again an exact Li₃ combination** (paper's tails method, `m3_tails_exact.py 5`; 14 Kummer/three-term/duplication/inversion relations; exact-vs-quadrature 0.0 at 60 dps, 2.5e-102 at 110 dps):
  W3⁽⁵⁾ = (−7/16+7√2/32)ζ(3) + (1+√2)Li₃(b) + (1−√2)Li₃(−b) − 2√2 Li₃(1/√2) + (−1/8−√2/24)λπ² − √2 λLi₂(b) + (−1/8+√2/12)Lπ² + 2L Li₂(b) + (−1/8+√2/24)λ³ + (√2/4)λ²L + (−1/2−√2/2)λL² + (5/6−√2/6)L³.

## 3. G-expression of W6 (`gexpr5.py` via `w4_symbol.log_product`; `gpl_w6.tsv`)

3158 terms, len(word)+k = 6 for all terms, alphabet {0, ±a, ±b}, base point 1, first letters only {0, −a} (no endpoint regularization needed). Evaluated with the acb Taylor evaluator (800 bits, cap 500): |G − quadrature| = **4.3e-101** (≥ 90-digit check), flint radius 1.7e-238.

## 4. 400-digit values (mpmath tanh-sinh on the v-form; two working precisions)

- W6: runs @520/@470 dps agree to **2.1e-421**; quad vs independent G-expression **2.6e-111**. `W6_400.txt`:
  W6 = −0.9076085816559868042386275171836577009010053662475073153459453747872290552306167366262240525183005019569583123474184260900325685558069756081943751142396145579672305516783069023591975635605782418141004363319860276506024120138276082834795123509635126427930308710073970036093582215646764314858632960459291891157972107854570495694258136789810634609641982780651239225446935354867669463684673053672424703075
- M5: runs @520/@470 dps agree to **2.5e-423**. `M5_400.txt`:
  M5 = 0.003161648710412257483507878511454280647014062547593368378732093712866724092945095130621400843393247057355686503030171199240441228324029373411335264512432749401798687387084452039624962371331452477259320689655410108805203773690659598986823320370987701071031535408464863416063743919282398530009541994013982636873731767709782203595599895036009485392030164706832039985196291937583640574106876904815460615758

## 5. The depth test at weight 6 (`w6_symbol.py`; full symbol in `symbol6.json`, 3158 G-terms, letters in span{X,Y})

Pairing outputs (verbatim):
```
control: 200 random shuffle products pair to zero: True
control: 100 Li6(S-unit) symbols pair to zero: True
Li_{4,2}(b,1) pairings (33, 42a,42b,42c): ['0', '0', '0', '-1/4']
Li_{3,3}(b,1) pairings (33, 42a,42b,42c): ['0', '0', '0', '3/4']
weight-5 control pairings ([X,P]xP, [Y,P]xP): ['84', '-38*sqrt2']  (expect 84 and -38 sqrt2)
W6 (3,3) pairing [X,P]x[Y,P]: 449/4*sqrt2
W6 (4,2) pairings ([X,[X,P]]xP, [Y,[X,P]]xP, [Y,[Y,P]]xP): ['-30', '15*sqrt2', '-9']
```
Independent recomputation from `symbol6.json` only (`planner6_pairing.py`, different code/RNG):
```
independent: 200 shuffle products pair to zero: True
independent: 100 Li6 symbols pair to zero: True
W6 (3,3) pairing: [0, 449/4]  (= 449 sqrt2/4)
W6 (4,2) pairings: [-30, 0], [0, 15], [-9, 0]  (= -30, 15 sqrt2, -9)
W5 control: [84,0], [0,-38]  (= 84, -38 sqrt2)
```
All mandatory controls pass: 200 shuffle products → 0; 100 Li₆(S-unit) symbols → 0; Li₄,₂(b,1) and Li₃,₃(b,1) → nonzero (−1/4 and 3/4 in the [Y,[Y,P]]∧P slot); the weight-5 run of the same code reproduces W5's (84, −38√2). Completeness of the test: verified exactly that rank(B3) = 2 = dim Lie³ and rank(B4) = 3 = dim Lie⁴ (Witt dims 2,1,2,3), so these four pairings are the *complete* (4,2)+(3,3) symbol-level cobracket. (As the brief notes, Λ³B2 = 0 since B2 = span{P}, so depth ≤ 2 is automatic at symbol level; the only question was whether depth two is *needed* — it is.)

## 6. The depth-two correction N6 (`w6_solve.py`; `n6_solution.json`)

Convention (verified numerically to 45 digits on test points): Li_{m,n}(x,y) = Σ_{p>q≥1} xᵖyᑫ/(pᵐqⁿ) = G₀(0^{m−1}, 1/x, 0^{n−1}, 1/(xy); 1).
Survey of 63 atoms (`atom_pairings.json`): every y = ±1 atom has zero (3,3) pairing (Li₃(±1) ∝ ζ(3) has zero symbol), so y ∉ {±1} is forced — hence y = 2 and y = 1/√2 below. Exact solve over Q(√2) (4×4 system; smallest of 2430 solvable 4-subsets):

**N6 = (−12 + 209√2)·Li₃,₃(b,1) + (−160 + 80√2)·Li₃,₃(1/√2,1) + (329√2/4)·Li₃,₃(b,2) + 120√2·Li₄,₂(b,1/√2)**

Direct exact verification (pairings recomputed from the summed atom symbols, not from the linear solve):
```
N6 (3,3): 449/4*sqrt2  target: 449/4*sqrt2  equal: True
N6 (4,2): ['-30', '15*sqrt2', '-9']  target: ['-30', '15*sqrt2', '-9']  equal: True
```
So W6 − N6 has zero (4,2) and (3,3) Lie cobracket at the symbol level. Numeric values (double series, cross-checked against the flint iterated-integral evaluator to ≤4.1e-73, radii ≤2.7e-121): N6 = 10.52854780888598960173545974996194338215864260058435519945…; W6 − N6 = −11.4361563905419764059740872671456010830596479668318625147959… (the classical part, closed form not required and not pursued). N6 is non-unique at symbol level (e.g. Li₃,₃(2,b) can replace Li₃,₃(b,2) with adjusted coefficients; the kernel is depth-two combinations that are symbol-classical).

## Evidence table (all residuals printed by code)

| check | residual | digits/precision |
|---|---|---|
| pieces(quadrature) sum − M5 | 1.0e-126 | 100-dps run |
| pieces(exact T, W3⁽⁵⁾) sum − M5 | 4.0e-100 | 110-dps run |
| Neville M5 − quadrature M5 | 2.2e-56 | independent method |
| T − (π²+6λ²)/24 | 0 (exact symbolic) + 0.0 | exact + 110 dps |
| W3⁽⁵⁾ exact − quadrature | 2.5e-102 | 110 dps |
| W6 G-expression − quadrature | 4.3e-101 (flint rad 1.7e-238) | ≥90-digit check |
| W6 two quadrature runs (520 vs 470 dps) | 2.1e-421 | 400+ digits saved |
| M5 two quadrature runs (520 vs 470 dps) | 2.5e-423 | 400+ digits saved |
| W6 quadrature − G-expression | 2.6e-111 | independent method |
| N6 atom values: series − integral evaluator | ≤ 4.1e-73 (rad ≤ 2.7e-121) | convention check |
| W5 control pairings of same code | exactly (84, −38√2) | exact |
| N6 pairings − W6 pairings | 0 | exact (QS2) |

## Producing files (all in `new_results/M5a_next_coefficient/`)

`ibp_reduce_k5.out`, `omega5.json` (Ω's; `omega5_parse.py` = auto-saved `step_19.py`), `quad5_100.json`, `m5_neville.json`, `m5_tails_exact.pkl`/`.out`, `gpl_w6.tsv`, `gpl_w6_check.json`, `hp5_results.json`, `W6_400.txt`, `M5_400.txt`, `symbol6.json`, `atom_pairings.json`, `n6_solution.json`, `n6_values.json`, `planner6_pairing.out`, `results_summary.json`, `progress.log`.
Scripts: `ibp_reduce.py`, `mk_integrands.py`, `qsqrt2.py` (copied from paper), `m3_tails_exact.py` (+`polylog_exact.py`, `w2_integrals.py`, `kummer.py`, `dilog_exact.py`, output paths redirected into this folder), `step_quad5.py`, `step_neville_m5.py`, `step_hp5.py`, `gexpr5.py`, `g_eval5.py`, `w6_symbol.py`, `w6_atoms.py`, `w6_solve.py`, `planner6_pairing.py`, `omega5_parse.py`; auto-saved verification snippets `step_42.py` (atom series/N6 value), `step_43.py` (exact-tails 110-dps check), `step_44.py` (Lie-rank checks), `step_46.py`/`step_47.py` (atom integral cross-check), `step_48.py` (summary).

## Rerun commands
```
cd "new_results/M5a_next_coefficient"
PY='<project>/.venv/bin/python'
nice $PY -B ibp_reduce.py 5            # Omega's -> ibp_reduce_k5.out
nice $PY -B omega5_parse.py            # parse -> omega5.json (exact reconstruction asserted)
nice $PY -B step_quad5.py              # pieces + M5 quadrature -> quad5_100.json
nice $PY -B step_neville_m5.py         # Neville check -> m5_neville.json
nice $PY -B m3_tails_exact.py 5        # exact T and W3^(5) -> m5_tails_exact.pkl/.out
nice $PY -B gexpr5.py                  # W6 G-expression -> gpl_w6.tsv
nice $PY -B g_eval5.py                 # G-expression vs quadrature -> gpl_w6_check.json
nice $PY -B step_hp5.py w6 520; nice $PY -B step_hp5.py w6 470
nice $PY -B step_hp5.py m5 520; nice $PY -B step_hp5.py m5 470   # -> hp5_results.json, W6_400/M5_400
nice $PY -B w6_symbol.py               # symbol + pairings + all controls -> symbol6.json
nice $PY -B w6_atoms.py                # atom pairing survey -> atom_pairings.json
nice $PY -B w6_solve.py                # N6 solve + exact verification -> n6_solution.json
nice $PY -B planner6_pairing.py        # independent pairing/control check
```
---

## Planner verification (W4 session, 2026-10-27)
- **Depth test, independent code** (`planner/planner_w6_pairing.py`, output `planner/planner_w6_pairing.out`): my own
  pairing code on the symbol of `gpl_w6.tsv`, built with the A3_W4_proof engine.
  - W6 gives (3,3) = 449√2/4 and (4,2) = (−30, 15√2, −9), identical to M5a's.
  - The controls all pass:
    - the weight-5 code reproduces (84, −38√2);
    - 200 shuffle products and 10 Li6 values give 0;
    - Li_{4,2}(b,1) and Li_{3,3}(b,1) give nonzero.
- **N6** (`planner/planner_check_n6.py`): the four atoms, with symbols built directly from their G-words, reproduce all
  four W6 pairings. So W6 − N6 has zero (3,3) and (4,2) parts, exactly.
- **M5 value, independent method:** K2_corner_decay/out/mk_pi_2.tsv (Cauchy/FFT at the square corner, dps 140) gives
  M5 = 0.00316164871041225748350787851145428064701406254759336837873209…, agreeing with M5a's 400-digit value in all 60
  digits shown.
- **Status:** M5's weight-6 part W6 needs a depth-two polylogarithm. This is VERIFIED at the symbol level (number level
  modulo the period conjecture), exactly as for M4's W5. T and W3⁽⁵⁾ are exact again (proved by the tails method).
  - Not done: closed forms of W6 − N6 (classical part) and W5″.
  - Typo in M5a's verdict line: "M6's top-weight part" should read "M5's".
