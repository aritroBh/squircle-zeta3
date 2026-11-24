# TOW4R — adversarial referee report on TOW4 (tower law to depth 8)

**Verdict: (a) tower law at depths δ ≤ 8 — PROOF COMPLETE (MINOR wording fixes). (b) weight-six law W₆(k) = Σₘ C(k−6,m−1)E_{5,m} for every k ≥ 6 — PROOF COMPLETE (the k = 6 level is an exact finite kernel check, as TOW4 states; the ±i pole vanishing is proved for k ≥ 7 by 5 zeros > degree 4 and at k = 6 exactly). (c) Task D reduction — correctly stated and verified at δ = 7,8; the all-depth law is honestly NOT claimed.** Every load-bearing certificate was re-executed by me in a byte-copied folder (`tow4_copy/`; TOW4's own folder never written). All residuals are EXACTLY 0 (exact rational arithmetic over Q(i) via sympy, or exact Fraction arithmetic in Q(√2,i); the only float anywhere is a 40-digit mpmath cross-check of τ(8), which is also proved symbolically). One displayed closed form is WRONG in both writeups (F(7,3)) — a typo, not a math error; the certificates and the actual family are correct.

## 1. Verdict table — every check I actually ran (all exact; residual = 0 means exactly zero)

| # | check (my rerun / my own code) | result | residual/evidence |
|---|---|---|---|
| 1 | `forced78.py` rerun: forced family δ≤8, closing F(7,4)=0, F(8,5)=0 | **PASS** | numerator 0 (39.5 s) |
| 2 | `prove_adj78_fast.py` all 9 adjacent identities (7,0..3),(8,0..4) | **PASS** | numerator 0 (2–82 s each) |
| 3 | `prove_E3b78.py` δ=7,8 (blindness) | **PASS** | numerator 0 (17.6/30.3 s) |
| 4 | `verify_tow4.py` full official runner | **10/10 PASS** | 413 s (log saved) |
| 5 | fast lcm-method vs TOW3 slow `sp.simplify` method on (7,3),(8,4),(8,3),(7,2) | **identical True** | both parities + both tops |
| 6 | vacuity: `check_zero` sanity + 3 perturbed feeds (sign flip; τ₁+1; fold constant i^{−6−δ}→i^{−4−δ}) | all **correctly FAIL** ("numerator == 0 : False", rc=1) | cannot pass vacuously |
| 7 | L cache == tow2lib.L(r) **symbolically**, r=1..8 (closes TOW3's r≤6 gap) | **equal** | diff 0 (r=8: 44 s) |
| 8 | τ₁(8)=−169/128+41i/128, τ₂ = conj, from the DL definition, exact Q(√2,i) + 40-digit mpmath; conj rules n=0..8 | **diff 0** | −1.3203125+0.3203125i ✓ |
| 9 | from-scratch L_r: genintegrand.S_poly + my own QSI series log-derivatives of A(v) at v=i (no tow2lib), 88 keys (d=6..18, incl. diagonal) × r=1..8 | **0 bad** | 1248 s |
| 10 | from-scratch a₇,a₈ (my Bell of my L) == T4.arho == towlib.rho (binteg initial data), 88 keys × 2 depths | **0 bad** | 973 s |
| 11 | RHO: a₇,a₈ vs TOWN1 towers at k=d, 717 cells each, **my own Q(√2,i) division** | **0 bad** | 3.9/22.6 s |
| 12 | `check_full78.py` rerun: full law δ≤8, 49953 cells + 222 diagonal-initial (k=2p) cells | **0 bad** | 21 s |
| 13 | empirical Newton N(7,s),N(8,s) from raw data only (off-diagonal, d=6..20, all s incl. top equations) | **846 cells, 0 bad** | exact Fractions |
| 14 | empirical adjacent N(δ,s) from raw data incl. fold (p=4..8; e00 ratios taken from data) | **45 cells, 0 bad** | |
| 15 | empirical E3 (diagonal blindness) from raw data, δ=7,8, p=4..12 | **18 cells, 0 bad** | |
| 16 | empirical diagonal ê-law from raw data, δ=7,8, p=3..13 | **99 cells, 0 bad** | |
| 17 | e00 ratios rA0=i(4p−3)/(p+1), rB0=(4p−1)(4p−3)/(p(p+1)), rB0d=2i(4p−1)/(p+1) vs data p=3..8; e00_{p,p−1} real, e00_pp anti-real | **all confirmed** | |
| 18 | **my own T1-in-k recursion** (fresh implementation of the TOW1 formula incl. fold) vs raw towers, degrees 12..20, k=d+1..52, δ=7,8 | **5628 cells, 0 bad** | (first run caught my own fold-scale bug: −i vs −p; fixed, rerun clean) |
| 19 | W6: k=6 vs quartic(7..11) **48/48** components; 5th differences vanish k=7..13 (48/48); ipos/ineg empty k=6..13; k=5 off the quartic at **15/48** | **all reproduced** | matches claims exactly |
| 20 | Task D telescoping identity F(δ,top+1)=(−1)^top[feed(U)−a_δ] with perturbed seed, symbolic δ=3,4,7,8; seed sensitivity | **identity 0; non-vacuous** | see §3 caveat |
| 21 | MASTER lemma feed(U)−a_δ==0 at δ=7,8 (inside verify_tow4 rerun) | **PASS** | |
| 22 | displayed closed forms vs pickles: a₇ full numerator, a₈ denominator + leading 256d¹²−10752d¹¹+(208384−3584i)d¹⁰, F(7,1),F(7,2),F(8,2),F(8,3) | **match** | **F(7,3) MISDISPLAYED** (§4.1) |

Digit counts: every residual above is an exactly-vanishing rational function/number. Data cells: ρ 2×717; law 49953+222; Newton 846; adjacent 45; E3 18; ê 99; T1-in-k 5628; kernels 3×48; from-scratch 88×8 (L) + 88×2 (a_δ).

## 2. Findings by the review items

**1. Genuinely symbolic — OK.** The depth-7,8 e-family is the forced ascending recursion from proved initial data a₇,a₈ (Bell over the L cache, which I proved equal to the proved tow2lib.L(r) symbolically for r=1..8) and τ to n=8 (proved). No fitting anywhere: the inherited `e_fits.pkl` (with the degenerate (7,3) fit) is never loaded by tow4lib or any proof; tow3lib.e enters only the CLOSE-sanity check (not load-bearing). All prove_*.py keep (d,m)/p symbolic; denominators are nonzero odd integers in the stated domains; the descending finite-induction uniqueness argument (CX8/CX9 form) is correctly invoked.

**2. ρ(7), ρ(8) — OK.** Re-derived from scratch: L_r = (log Φ)^{(r)}(i)/r! with Φ(v) = v^{2d−1}S_m(t(v))(v+i)^{2−2d}/((v²−a²)(v²−b²)), using genintegrand.S_poly, my own exact QSI series arithmetic (log via f′/f), and my own check t′(i)=0, t^{(r)}(i)/r! = i^{r+1}/4 (= Lemma T). My L_r == tow2lib.L(r) (88 keys × r=1..8); my Bell == T4.arho == binteg ratios (88 keys × 2); 717 cells each vs TOWN1 with my own division. R₅ correctly not needed (partitions of r ≤ 9 into parts ≥ 2 have ≤ 4 parts); the "ρ(10) needs R₅" boundary is correct.

**3. Newton identities, fold, ê, E3 — OK.** Non-adjacent: forced recursion correct; closing F(7,4)=0, F(8,5)=0 exact (rerun). Adjacent: 9/9 exact; the new fast method computes DEN·(lhs−Σterms) exactly (lcm of monic Q(i)-irreducible denominator factors; `assert r==0` on exact division guards it), validated against the slow method on 4 cases (both parities, both tops), and cannot pass vacuously (3 perturbations fail). The feed structure, ê alternating conversion, fold constants (−i)^{−4−t}i^{−6−δ} (i^{4p}=1 for integer p — I verified the reduction against towlib's i^{−mu}i^{mm−2}), rA0/rB0/rB0d, and Dc = −conj(Δ) (anti-real e00_pp) are byte-identical in structure to the twice-refereed TOW3 scripts (diffs shown) and are all confirmed empirically at the new depths (rows 13–18).

**4. Logic chain to (b) — OK.** Degree-6 towers to depth 8 = 2·6−4 hold for k ≥ 7 (off-diagonal k ≥ 6; diagonal (3,3) stable k ≥ 7). WTS Lemma 3 (statement read and checked) ⇒ degree-5 kernel residues polynomial of degree ≤ 4 for k ≥ 7; Lemma 2 (fixed real-linear image, conjugates included — CX8 #10 phrasing carried forward) ⇒ W₆(k) quartic for k ≥ 7. ±i simple-pole coefficients are residues ⇒ quartics, vanish at k=7..11 (5 = d points) ⇒ identically for k ≥ 7; k=6 by the exact 48/48 finite check (±i lists empty at k=6). k=5 confirmed exceptional (15/48). Weights polynomial in k for ALL k ≥ w: w=2,3,4,5 (previous), **w=6 (here)**; w≥7 open (needs δ≤10). Top-coefficient table: zigzag tops are now the PROVED actual coefficients for δ ≤ 8 (A₂..A₅ = 1,2,5,16; adjacent tops (7,3),(8,4) newly proved); for δ ≥ 10 the top equations are proved but numerators remain value-forced/data-verified. Paper `prop:w6k` (6 ≤ k ≤ 52, "under review" comment) is accurate and may become "every k ≥ 6".

**5. Independent data — OK.** Rows 9–11, 13–18 use only raw TOWN1 towers, WTS binteg, and genintegrand — no tow2lib/tow4lib formulas — at fresh degrees/levels beyond TOW2/TOW3's checks.

**6. Task D — OK as stated.** (1)⟺(2) is telescoping (verified symbolically with a perturbed seed, hence non-vacuous); (2)⟺(3) is the definition of the formal x=−1 specialization of (1+x)Û = ρ + xΛ(Û); (1)⟺(4) is the same telescoping read descending against the proved tops. Certificate feed(U)−a_δ ≡ 0 at δ=7,8 reruns green. The obstruction (U(δ) self-reference at dp=δ; no all-depth proof) is honest — see §3.

## 3. Caveat (sharpening, not a gap): even-depth seed blindness

Perturbing the seed a_δ → a_δ+ε (lower family fixed): d(F(δ,top+1))/dε is **nonzero at odd δ = 3,5,7** (δ=7: 105/(16d⁴−160d³+560d²−800d+384)) but **exactly 0 at even δ = 4,6,8** (verified symbolically). So at even depth the closing identity — equivalently the master x=−1 identity — is blind to a_δ itself: the U(δ) self-reference in the feed exactly cancels −a_δ, leaving a constraint purely on the lower-depth family. (i) The δ=8 closing check alone would NOT have caught a wrong a₈; ρ(8) rests entirely on Task A's initial-data derivation — which is proved (R₁..R₄ → L → Bell) and which I verified from scratch (rows 9–11), so there is no gap. (ii) For the all-depth problem this mildly weakens TOW4's "self-referential" obstruction at even δ (the seed persists only at odd δ). Worth one sentence in proof_depth8.md §5.

## 4. Wording fixes (MINOR; none load-bearing)

1. **F(7,3) is misdisplayed in BOTH result.md and proof_depth8.md** as "F(7,3) = −2i(d−5)/Π₃ (zigzag, A₄ = 5)". The proved value (forced_F78.pkl; consistent with A₄=5 and with the F(top)==zigzag sanity, which I reproduced) is **F(7,3) = +5i(d−5)/Π₃**. As written the display contradicts the zigzag formula and the (7,3) adjacent certificate. (The stale e_fits (7,3) fit is degenerate and matches neither.)
2. **Task D claim (3) parenthetical** "x⁰ = value identities, x^{s≥1} = slopes, x^{top+1} = tops" is off by one vs CX8 #7 for (1+x)Û = ρ + xΛ(Û): x⁰ is the initial condition e(δ,0)=ρ(δ); x¹ is the value identity; x^s is N(δ,s−1); x^{top+1} is the top. The x=−1 specialization itself is stated correctly.
3. **Domains:** "adjacent needs 4p−6 ≥ δ (p ≥ 3)" is the δ≤6 specialization; at δ = 7,8 adjacent keys and E3 need **p ≥ 4** (diagonal tower needs 4p−4 ≥ δ). State the per-depth bound.
4. **Stale docstring:** TOW4's byte-copied tow2lib.py still says "R_3, R_4 fitted+verified (prove_S34.py pending)" though prove_S34.py supplies the proof (CX8 #8 recurs; proof_depth8.md states it correctly).
5. `verify_tow4.py` greps subprocess stdout without checking return codes (LAW/ADJ/E3) — robustness nit only; my reruns checked return codes.
6. Mark inherited `e_fits.pkl` (7,·),(8,·) as unused candidates (TOW3R fix #5 pattern).
7. Add the even-depth seed-blindness sentence (§3) to Task D.

## 5. Deliverables in `new_results/TOW4R_review/`

`result.md` (verdict-first), `verify_tow4r.py` (10 items, < 10 min, PASS/FAIL: TAU8; L r=1..8; forced78 closing; ADJ tops+(8,3); E3 7,8; non-vacuity; LAW 49953+222; RHO 2×717 own division; W6; MASTER telescoping), `verify_tow4_rerun.log` (10/10, 413 s), `progress.log`, `tow4_copy/` (path-patched copies), `tow2_orig/`, `tow3_orig/` (diff sources). Step files that produced the results (harness-written, same folder): step_01–24 (setup; reads; CX8/CX9/TOW2R/TOW3R/WTS/paper analysis; diffs prove_E3b56→78 and prove_adj56→78_fast; F(7,3) misdisplay; a₇/a₈ display checks; τ(8) exact+40-digit; L r=1..8 symbolic), step_25–44 (from-scratch L_r 88 keys; a₇/a₈ vs arho/binteg; forced78 rerun; 9 ADJ + 2 E3 reruns; vacuity tests; fast-vs-slow cross-validation), step_45–57 (empirical Newton 846; adjacent 45; E3 18; ê 99; e00 ratios; RHO 2×717 own division; check_full78 rerun; W6 kernels; telescoping + seed sensitivity δ=3..8), step_58–70 (own T1-in-k 5628 cells incl. my fold-scale bug fix; verify_tow4 full rerun 10/10; e_fits provenance).

**Bottom line:** the tower law at δ ≤ 8 (all keys, all degrees, all k; off-diagonal k ≥ d; diagonal k ≥ 2p+1 with k = 2p exceptional) and W₆(k) = Σ_{m=1..5} C(k−6,m−1) E_{5,m} for every k ≥ 6 stand as proved. No mathematical errors found; the fixes above are presentational.