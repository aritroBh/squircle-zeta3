
result = r'''# TOW4R — adversarial referee report on TOW4 (tower law to depth 8)

**Verdict: (a) tower law at depths δ ≤ 8 — PROOF COMPLETE (MINOR wording fixes). (b) weight-six law W₆(k) = Σₘ C(k−6,m−1)E_{5,m} for every k ≥ 6 — PROOF COMPLETE (the k=6 level is an exact finite kernel check, as stated). (c) Task D reduction — correctly stated and verified at δ=7,8; the all-depth law is (honestly) not claimed.** Every load-bearing certificate was re-executed by me in a byte-copied folder (`tow4_copy/`); all residuals are EXACTLY 0 (exact rational arithmetic over Q(i) / exact Fraction Q(√2,i); the only floating point anywhere is a 40-digit mpmath cross-check of τ(8), which is also proved symbolically). One displayed closed form is WRONG in both writeups (F(7,3); the proved value is +5i(d−5)/Π₃) — a typo, not a math error; the certificates and the actual family are correct.

## 1. What I ran (all exact; residual = 0 means exactly zero)

| # | check (my rerun / my own code) | result |
|---|---|---|
| 1 | `forced78.py` rerun: forced family δ≤8, closing F(7,4)=0, F(8,5)=0, numerator 0 | **PASS** (39.5 s) |
| 2 | `prove_adj78_fast.py` all 9 adjacent identities (7,0..3),(8,0..4), numerator 0 | **PASS** (45/19/6/2/82/37/14/4/1 s) |
| 3 | `prove_E3b78.py` δ=7,8, numerator 0 | **PASS** (17.6/30.3 s) |
| 4 | `verify_tow4.py` full official runner | **10/10 PASS** (413 s) |
| 5 | fast lcm-method vs TOW3 slow `sp.simplify` method: (7,3),(8,4),(8,3),(7,2) | **identical True** (3/1/17/43 s) |
| 6 | vacuity: `check_zero` sanity + 3 perturbed feeds (sign flip; τ+1; fold constant) | all **correctly FAIL** |
| 7 | L cache == tow2lib.L(r) symbolically, r=1..8 (the r=7,8 gap left at TOW3) | **equal** (r=8: 44 s here) |
| 8 | τ1(8),τ2(8) from the DL definition, exact in Q(√2,i), + 40-digit mpmath; conj rules n=0..8 | **diff 0** |
| 9 | from-scratch L_r: genintegrand.S_poly + my own exact QSI series implementation of log A derivatives at v=i (no tow2lib), 88 keys × r=1..8 | **0 bad** (1248 s) |
| 10 | from-scratch a₇,a₈ (my Bell of my L) == T4.arho == towlib.rho (binteg initial data), 88 keys × 2 depths | **0 bad** (973 s) |
| 11 | RHO: a₇,a₈ vs TOWN1 towers at k=d, 717 cells each, with my own Q(√2,i) division | **0 bad** (3.9/22.6 s) |
| 12 | `check_full78.py` rerun: full law δ≤8 vs data, 49953 cells + 222 diagonal-initial cells | **0 bad** (21 s) |
| 13 | empirical Newton N(7,s),N(8,s) from raw data only (off-diagonal, d=6..20, all s incl. tops) | **846 cells, 0 bad** |
| 14 | empirical adjacent N(δ,s) from raw data incl. fold (p=4..8, e00 ratios taken from data) | **45 cells, 0 bad** |
| 15 | empirical E3 (diagonal blindness) from raw data, δ=7,8, p=4..12 | **18 cells, 0 bad** |
| 16 | empirical diagonal ê-law from raw data, δ=7,8, p=3..13 | **99 cells, 0 bad** |
| 17 | e00 ratios rA0=i(4p−3)/(p+1), rB0=(4p−1)(4p−3)/(p(p+1)), rB0d=2i(4p−1)/(p+1) vs data p=3..8; e00_{p,p−1} real, e00_pp anti-real | **all confirmed** |
| 18 | my own T1-in-k recursion (fresh implementation of the TOW1 formula incl. fold) vs raw towers, degrees 12..20, k=d+1..52, δ=7,8 | **5628 cells, 0 bad** |
| 19 | W6: k=6 vs quartic(7..11) 48/48 components; 5th differences vanish k=7..13 (48/48); ipos/ineg empty k=6..13; k=5 off the quartic at 15/48 | **all reproduced** |
| 20 | Task D telescoping identity F(δ,top+1)=(−1)^top[feed(U)−a_δ] with perturbed seed, symbolic (δ=3,4,7,8); seed sensitivity d/dε | **identity 0; non-vacuous** |
| 21 | MASTER lemma (feed(U)−a_δ==0 at δ=7,8) — rerun inside verify_tow4 | **PASS** |
| 22 | displayed closed forms vs forced_F78.pkl/a7_closed.pkl/a8_closed.pkl: a₇ full numerator, a₈ denom+leading, F(7,1),F(7,2),F(8,2),F(8,3) | **match**; **F(7,3) MISDISPLAYED** (below) |

Digit counts: every residual above is an exactly-vanishing rational function/number (sympy Q(i) numerator 0, or exact Fraction Q(√2,i) difference 0). Data cells: ρ 2×717; law 49953+222; Newton 846; adjacent 45; E3 18; ê 99; T1-in-k 5628; kernels 3×48; from-scratch 88×8 (L) + 88×2 (a_δ).

## 2. Findings by the review items

**1. Genuinely symbolic — OK.** The e-family at depths 7,8 is the forced ascending recursion from the proved initial data a₇,a₈ (Bell assembly over the L cache, which I proved equal to the proved tow2lib.L(r) symbolically for r=1..8) and τ extended to n=8 (proved). No fitting anywhere: `e_fits.pkl` (with its degenerate (7,3) fit) is inherited from TOW3 and is NOT used by any TOW4 proof (tow4lib never loads it). tow3lib.e enters only the CLOSE-sanity check (not load-bearing). All prove_*.py treat (d,m) resp. p as free symbols; denominators are nonzero odd integers in the stated domains. The descending finite-induction uniqueness argument (CX8/CX9 form) is correctly cited and applies.

**2. ρ(7), ρ(8) — OK.** I re-derived L_r = (log Φ)^{(r)}(i)/r! from the analytic factor Φ(v) = v^{2d−1}S_m(t(v))(v+i)^{2−2d}/((v²−a²)(v²−b²)) from scratch (genintegrand S_poly, my own QSI series arithmetic, my own Faa di Bruno via series f′/f; t′(i)=0, t^{(r)}(i)/r! = i^{r+1}/4 matches Lemma T): my L_r == tow2lib.L(r) at 88 keys × r=1..8; my Bell assembly == T4.arho == binteg initial-data ratios at 88 keys × 2 depths; and 717 cells each vs TOWN1 at k=d with my own division. R₅ is correctly NOT needed (partitions of r≤9 into parts ≥2 have ≤4 parts); the honest boundary (ρ(10) needs R₅) is correct.

**3. Newton identities, fold, ê, E3 — OK.** Non-adjacent: the forced recursion is correct and the closing conditions F(7,4)=0, F(8,5)=0 hold as exact rational identities (rerun). Adjacent: the 9 identities hold; the new fast method is a sound exact computation of DEN·(lhs−Σterms) (I analyzed the lcm/monic/content logic line by line; the `assert r==0` on exact division guards it), is validated against the slow method on 4 cases (both parities + both tops), and CANNOT pass vacuously (3 perturbed feeds all fail). The feed formulas, ê conversion Σ_t(−1)^t ê(s+t), fold constants (−i)^{−4−t}i^{−6−δ} (= i^{4p−4−t→} reduced mod 4, valid for integer p), rA0/rB0/rB0d, and Dc = −conj(Δ) (anti-real e00_pp) are byte-identical in structure to the TOW3 scripts both referees validated, and I confirmed all of them empirically against raw data at the NEW depths (items 13–17 above), plus the T1 recursion itself at 5628 fresh cells (item 18).

**4. Logic chain to (b) — OK.** Degree-6 towers to depth 8 = 2·6−4 are proved for k ≥ 7 (off-diagonal k ≥ d = 6 ⊇ k≥7; diagonal (3,3) stable k ≥ 2p+1 = 7). WTS Lemma 3 (statement checked) with d=5, n=6 gives degree-5 kernel residues polynomial in k of degree ≤ 4 for k ≥ 7; WTS Lemma 2 (fixed real-linear image; conjugates included — CX8 #10 phrasing correctly carried forward) gives W₆(k) quartic for k ≥ 7. The ±i simple-pole coefficients are residues, hence quartics, vanish at k=7..11 (5 points = d), hence identically for k ≥ 7. k=6 (the exceptional level of (3,3), outside Lemma 3) is settled by the exact finite check: 48/48 components on the quartic, ±i lists empty at k=6. k=5 confirmed exceptional (15/48 mismatch). Weights polynomial in k for ALL k ≥ w: w=2,3,4,5 (previous), **w=6 (here)**; w≥7 open (needs δ≤10). Top-coefficient table (LEAD1/LEAD2): the zigzag tops are now the PROVED actual coefficients for δ ≤ 8 (A₂..A₅ = 1,2,5,16; adjacent tops (7,3),(8,4) newly proved); for δ ≥ 10 the top equations remain proved but numerators are only value-forced/data-verified. Paper `prop:w6k` currently says 6 ≤ k ≤ 52 with an explicit "under review" comment — accurate; it may become "every k ≥ 6".

**5. Independent data — OK.** Items 9–11, 13–18 are all independent of TOW4's symbolic machinery (raw TOWN1 towers, WTS binteg, genintegrand), at fresh degrees/levels beyond TOW2/TOW3's checks.

**6. Task D — OK as stated.** The reduction (1)⟺(2) is the telescoping of the recursion (I verified the identity F(δ,top+1) = (−1)^top[feed(U)−a_δ] symbolically with a perturbed seed, so it is not vacuous); (2)⟺(3) is the definition of the formal x=−1 specialization of (1+x)Û = ρ + xΛ(Û); (1)⟺(4) is the same telescoping read descending (V(δ) against the proved tops). The certificate feed(U)−a_δ==0 at δ=7,8 reruns green. The obstruction (self-reference of U(δ) at dp=δ; no all-depth proof) is honestly stated.

## 3. Caveat (sharpening, not a gap): even-depth seed blindness

With the seed perturbed, a_δ → a_δ+ε (lower family unperturbed): d(F(δ,top+1))/dε is **nonzero at odd δ (3,5,7)** but **exactly 0 at even δ (4,6,8)** (verified symbolically; e.g. δ=7: 105/(16d⁴−160d³+560d²−800d+384); δ=8: 0). So at even depth the closing identity — equivalently the master x=−1 identity (2) — is blind to the initial data a_δ itself: it is a constraint purely on the lower-depth family (the U(δ) self-reference in the feed exactly cancels −a_δ). Consequences: (i) the δ=8 closing check alone would NOT have caught a wrong a₈ — ρ(8) rests entirely on Task A's initial-data derivation, which is proved (R₁..R₄ → L → Bell) and which I verified from scratch (items 9–11), so there is no gap; (ii) for the all-depth problem, at even δ the identity (2) effectively involves only the lower family, which mildly weakens TOW4's stated "self-referential" obstruction at even depths (at odd δ the seed persists). Worth one sentence in proof_depth8.md §5.

## 4. Wording fixes (MINOR; none load-bearing)

1. **F(7,3) is misdisplayed in BOTH result.md (Task B) and proof_depth8.md (§2)** as "F(7,3) = −2i(d−5)/Π₃ (zigzag, A₄ = 5)". The proved value (forced_F78.pkl; consistent with A₄=5 and with the F(top)==zigzag sanity, which I reproduced) is **F(7,3) = +5i(d−5)/Π₃**. As written the display contradicts both the zigzag formula and the adjacent (7,3) certificate.
2. **Task D, claim (3) parenthetical:** "x⁰ = value identities, x^{s≥1} = slopes, x^{top+1} = tops" is off by one relative to CX8 #7 for (1+x)Û = ρ + xΛ(Û): the x⁰ coefficient is the initial condition e(δ,0) = ρ(δ); x¹ is the value identity N(δ,0); x^s is N(δ,s−1); x^{top+1} is the top. The x=−1 specialization itself is stated correctly.
3. **Domains:** "adjacent needs 4p−6 ≥ δ (p ≥ 3)" is the δ≤6 specialization; at δ = 7,8 the adjacent keys and E3 need **p ≥ 4** (and the diagonal's own tower needs 4p−4 ≥ δ). State the per-depth bound.
4. **Stale docstring:** TOW4's byte-copied tow2lib.py still says "R_3, R_4 fitted+verified (prove_S34.py pending)" although prove_S34.py supplies the proof (CX8 #8 recurs). proof_depth8.md states it correctly.
5. `verify_tow4.py` greps subprocess stdout without checking return codes (LAW/ADJ/E3 lemmas) — robustness nit only; I reran everything with returncode checks.
6. Mark the inherited `e_fits.pkl` (7,·),(8,·) entries as unused candidates (TOW3R fix #5 pattern); the (7,3) fit is degenerate (b=0) and agrees with neither the misdisplayed nor the true value.
7. Add the even-depth seed-blindness sentence (§3 above) to Task D.

## 5. Deliverables in `research/08_theory/TOW4R_review/`

`result.md` (this file), `verify_tow4r.py` (10 items, <10 min, PASS/FAIL), `verify_tow4_rerun.log` (10/10), `progress.log`, `tow4_copy/` (byte copies + path-patched runners; TOW4's own folder never written), `tow2_orig/`, `tow3_orig/` (for diffs). Harness step files that produced the results: step_01–24 (setup, reads, diffs, F(7,3) analysis, a₇/a₈ display checks, τ(8), L r=1..8), step_25–44 (from-scratch L_r 88 keys; a₇/a₈ vs arho/binteg; E3/adjacent structural diffs; forced78 rerun; 9 ADJ + 2 E3 reruns; vacuity tests; fast-vs-slow), step_45–57 (empirical Newton 846; adjacent 45; E3 18; ê 99; e00 ratios; RHO 2×717 own division; check_full78 rerun; W6 kernels; telescoping + seed sensitivity δ=3..8), step_58–70 (own T1-in-k 5628 cells incl. fold-scale correction; verify_tow4 full rerun 10/10; e_fits provenance).

**Bottom line:** the tower law at δ ≤ 8 (all keys, all degrees, all k; off-diagonal k ≥ d, diagonal k ≥ 2p+1 with k=2p exceptional) and W₆(k) = Σ_{m=1..5} C(k−6,m−1)E_{5,m} for every k ≥ 6 stand as proved. No mathematical errors found; the fixes above are presentational.'''
import os
work = "<project>/research/08_theory/TOW4R_review"
open(work + "/result.md", "w").write(result)
print("result.md written", len(result.splitlines()), "lines")
