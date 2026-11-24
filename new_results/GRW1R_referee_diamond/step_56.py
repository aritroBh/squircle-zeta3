result = r'''# GRW1R — Adversarial referee report on GRW1 (diamond end: growth law of Taylor coefficients, Borel summability)

Referee: independent re-derivation and arb/mpmath verification. All numbers below were printed by the
referee's own code (step files step_*.py, master script verify_grw1r.py, 40.7 s). GRW1's labels:
PROVED = by hand; CERTIFIED = arb; NUMERICAL = mpmath.

## 1. Verdict table

| Claim | Statement (short) | Verdict |
|---|---|---|
| A | Laplace form I(η) = ∫₀^∞ e^{−u} H(ηu) K(η,u) du, η > −1 | **PROVED** |
| B | f_n = Σ_{j=0}^n H_{n−j} M_{n,j} exactly; M_{n,0} = n! η_D(n); |M_{n,j}| ≤ (n−j)!/(j!(j+1)^{n−j+1}) | **PROVED** |
| C | Darboux transfer for H: E* = 5/24 on the slits; |R_n| ≤ 4.45 ρ^{−n} n^{−7/2} + 135 (2.9ρ)^{−n} (n ≥ 10) | **PROVED AFTER REPAIRS** (E* now certified; phase typo fixed; remainder derivation defective — repaired constant 7.06; stated 4.45-bound is true, verified n ≤ 240) |
| D | Uniformity: |f_n − H_n n! η_D(n)| ≤ C_H e^{2ρ} ρ^{−n} n! 2^{−n−1}/n | **PROVED** (given C as repaired) |
| E | Growth law, square: f_n = √2 cos(π(3−2n)/4) n! (2/π)^n n^{−3/2} (1 + (3−3π(−1)^n)/(8n) + O(n^{−2})); rhombus law with θ₀; radius zero | **PROVED AFTER REPAIRS** (inherits C's repairs; the O(n^{−2}) term *is* proved) |
| F | Borel summability along arg ζ = 0 and π; S(η) = I(η) for every η > −1 | **PROVED AFTER REPAIRS** (Lemma 6's J₀-series "strip extension" is false — repaired via Lemma 6b; 6b's contour arithmetic to be spelled out; two cosmetic fixes) |
| G | Singular directions; decomposition of B in |Im ζ| < π | **PROVED** (with GRW1's own caveat: possible R₂ line-boundaries on |Im ζ| = (j+1)π/2, j ≥ 1, not pinned down; "exactly ±π/2 singular directions" holds modulo that caveat; arg ζ = 0, π are proved nonsingular) |
| H | Stokes jump 2√π e^{−3iπ/4} ε^{1/2} e^{−π/(2ε)} (1 − 3iε/4 + O(ε²)) | **PROVED** (given the proved singularity structure); asymptotics verified, next coefficient −25/32 extracted |

## 2. Line-by-line ledger of proof.md

**Lemma 1 (A). CORRECT.** The chain P = 4r(1+η)∫₀¹[t(1−t)]^{η/2}g(ηl)dt → fold at l = 0 → the displayed
integral is right; convergence for η > −1 fine. Referee: independent check at η = 0.37, Laplace form vs
perimeter integral, |diff| = 7.8e-62 (step_28); planner had 90 digits at 4 points. Lemma 1(a), 1(b):
CORRECT (|k_j| ≤ e^{−ju}/j! since 0 < log(1+e^{−u}) ≤ e^{−u}; note k_j = (−1)^j|k_j| — sign is as used).
Lemma 1(c): CORRECT (|H(z)| ≤ √2(1+e^{−2Re z})^{1/2} ≤ 2 max(1, e^{−Re z}) in the square).

**Lemma 2 (B). CORRECT.** Exact term-by-term identity; the dominator for ∂_ηⁿ|₀∫ = ∫∂_ηⁿ|₀ is correctly
built (Cauchy on disks of radius π/4 avoiding the branch points at ±iπ/2). Referee: own quadrature,
worst relative error 8.7e-62 at n = 1..5, plus n ≤ 55 to quad tolerance; M(1,1) = (ln2−1)/2 exact;
M_{n,0} = n!η_D(n) and the M-bound spot-verified (worst ratio 1.0 at j = 1, bound tight). (step_22–24)

**Lemma 3.1. CORRECT** (constants c ∈ {A√2/r, B√2/r}; the local form w^{1/2}(1 − w/2 + O(w²))).

**Lemma 3.2 (C). GAP → repaired.** Three findings:
1. *"E* = 5/24 certified on a 4000-point grid" is not a proof.* REPAIRED: E* = 5/24 is now a theorem.
   Referee proof (step_03, step_08): ψ(s) := |E(is)|² = |e^{−is/2}√(sin s/s) − 1 + is/2|²/s⁴ is a power
   series Ψ(s²) = Σc_m s^{2m} with c_0 = (5/24)², c_1 = −25/13824 < 0, c_2 = 11147/232243200 (exact
   rationals). On [0, 0.1]: Ψ(u) ≤ (5/24)² + (c_1 + (c_2+T)u)u with T = M_E(3)²/(9³(1−u/9)),
   M_E(3) ≤ 1.1900048 (rigorous crude bound |E| ≤ (√((1+e⁶)/6)+1+1.5)/9 on |w| ≤ 3 < π), and
   c_1 + (c_2+T)·0.01 ≤ −1.7885e-3 < 0 (arb) ⇒ ψ ≤ (5/24)². On [0.1, 1.5]: the same series (45 exact
   terms + Cauchy tail 1e-24) evaluated on arb intervals (h = 1e-3, 700 intervals): max ψ ≤ 0.0433847019.
   On [1.5, 1.9ρ]: direct arb interval evaluation of the closed form (744 intervals): max ψ ≤ 0.0399729017.
   All < (5/24)² = 0.0434027778. **Sup |E| on the slit = 5/24, attained only at s = 0. CERTIFIED.**
2. *Phase typo in the printed Hankel display.* GRW1 prints
   ℋ₊ + ℋ₋ = (2/π)ρ^{1/2−n}B(3/2,n−1/2)·2Re[e^{iπ/4}i^{−n}] − (1/π)ρ^{3/2−n}B(5/2,n−3/2)·2Re[i e^{iπ/4}i^{−n}].
   With these literal phases the main term has the WRONG sign pattern at even n ((+,+,−,−) instead of
   (−,+,+,−)) and the second term the wrong phases at odd n. The correct phases (referee derivation,
   orientation fixed so that H_n = A_n + B_n + (E-terms) + (tails) + (circle) holds with |H_n − A_n − B_n|
   ≤ E-term bound at every n tested, odd and even, n = 11..240):
     A_n = −(2/π)·2cos(nπ/2 + π/4)·ρ^{1/2−n}·B(3/2, n−1/2),
     B_n =  (2/π)·cos(nπ/2 − π/4)·ρ^{3/2−n}·B(5/2, n−3/2).
   (Equivalently, 2Re[e^{iπ/4}i^{−n}] must read 2Re[−e^{iπ/4}i^n] etc.; GRW1's final two-term law is
   exactly what these phases give — so this is a writeup typo, not a math error.) MINOR FIX.
3. *The remainder derivation is defective.* The printed intermediate
   "Γ(n−5/2)/Γ(n+1) ≤ 1.63·n^{−7/2} for n ≥ 10" is **false at n = 10** (true ratio 1.6306839; the
   resulting E-term bound 4.45·u survives anyway: certified sup 4.4454 at n = 10). The clause
   "the O(n^{−2}) truncations of the two Beta expansions are absorbed in the same 4.45" is **wrong as
   accounting**: |A_n + B_n − law2(n)| is a *separate, additive* piece, sup 2.7882·u at n = 11
   (u = ρ^{−n}n^{−7/2}), → 2.02·u as n → ∞. Honest assembly (arb-certified at each integer 10 ≤ n ≤ 500,
   sup at n = 11; rigorous Stirling tail with remainder bounds for n ≥ 500 giving ≤ 5.19):
     |R_n| ≤ |A+B−law2| + (4/π)E*ρ^{5/2−n}B(7/2,n−5/2) + tails + 135(2.9ρ)^{−n}
          ≤ **7.06**·ρ^{−n}n^{−7/2} + 135(2.9ρ)^{−n},  all n ≥ 10.   (certified sup 7.0523 at n = 11)
   The *stated* bound 4.45·u + 135(2.9ρ)^{−n} is **true** (exact rational H_n: worst ratio 0.9864 at
   n = 18 on 10 ≤ n ≤ 240; the 135(2.9ρ)^{−n} term rescues n = 10–17; asymptotically |R_n|/u → ~3.98 < 4.45),
   but it is not proved by the printed argument; proving 4.45 would need the E-integral expanded one
   term further (its leading part is explicit), not done here. With 7.06 every downstream constant
   survives: |H_n| ≤ (1 + 0.00817 + 7.06/190²)·|lawH| ≤ 1.032·lawH and ≥ 0.968·lawH for n ≥ 190 ✓.

**Lemma 3.3. CORRECT** (given repaired 3.2). Verified: max |H_n|(π/2)^n n^{3/2} = 5.8739238 at n = 3
(exact rationals, n ≤ 220); min |H_n/lawH| = 0.87454216 at n = 10 on 10..200; two-term residual
|H_n − law2|·ρⁿn^{7/2} ∈ [3.98, 4.58] for even n ≤ 240 (odd n smaller, ≈ 0.28–0.66; sign oscillation).
The sign check −cos((2n+1)π/4) = pattern (−,+,+,−) ✓ equals the paper's cos(π(3−2n)/4) pattern.
Note: the pattern is a property of the law; actual signs of f_n match it for 6 ≤ n ≤ 59 measured
(deviations only at n = 0, 1, 4, 5 — non-asymptotic) and provably for all n ≥ 10.

**Lemma 3.4. CORRECT** for the leading term and the limsup (bracket = r cos(ψ_n ± φ*) by parity; the
nonvanishing-orbit argument is sound: orbit of step 2θ₀ ∈ (0, π), 2θ₀ ≠ π, has ≥ 3 points mod 2π,
can't lie in the 2-point zero set of the cosine). The per-case explicit E*_j constants are NOT certified
(only numerical, 3 cases) — but the paper's thm:zero needs only the limsup, which uses soft bounds.
Referee spot-check of the bracket for A=1,B=2, A=1,B=3, A=2,B=3: max |bracket|/r ≈ 1 ✓. (step_30)

**Lemma 4 (D). CORRECT.** Every inequality verified: |M_{n,j}| ≤ (n−j)!/(j!(j+1)^{n−j+1}) (tight, worst
ratio 1.0 at j = 1); |H_m| ≤ C_Hρ^{−m}m^{−3/2} (m ≥ 1) and |H_0| = 2 ≤ C_H for the m = 0 term (this
convention should be made explicit — MINOR writeup note); (n−j)! ≤ n!/n; (j+1)^{n−j+1} ≥ 2^{n−j+1};
Σ_{j≥1} ρ^j/(j! 2^{n−j+1}) ≤ 2^{−n−1}e^{2ρ} (slack only 1.045 at n = 120). Referee numerics with the
450-digit f_n (step_12, step_31):
  n=20: |tail| = 2.8433e5 vs bound 9.427e8; n=40: 1.7789e24 vs 1.802e28; n=80: 3.450e74 vs 1.028e79;
  n=120: 1.502e134 vs 8.327e138. Relative to main: 9.2e-8, 3.9e-14, 1.7e-26, 1.0e-38 — the tail is
  2^{−n}-exponentially smaller, **the uniformity step the paper left unwritten is correct and proved.**

**Theorem 1 (E). PROVED AFTER REPAIRS (inherits C).** Square: the two-term law is right and the
O(n^{−2}) is really proved: |R_n|/|lawH| ≤ (7.06 + 135·2.9^{−n}n^{7/2})·n^{−2} ≤ 17.2·n^{−2} for
n ≥ 10 (sup of the circle term 10.15 at n = 10), and the K-tail is O(n^{1/2}2^{−n}) relative (Lemma 4),
η_D(n) = 1 + O(2^{−n}) — both absorbed. Per-class Richardson on the 450-digit f_n (n ~ 152–176):
constant 0.9999977–0.9999978 (all four classes), c_1 = −0.80202 (pred (3−3π)/8 = −0.80310) even classes,
+1.55414 (pred (3+3π)/8 = +1.55310) odd classes; two-term residual ×n² → −4.05 (even), +0.43 (odd).
Radius zero (thm:zero): PROVED (square: |cos((2n+1)π/4)| = √2/2 never vanishes ⇒ |f_n/n!|^{1/n} → 2/π
directly; rhombus: limsup = 1/ρ via Lemma 3.4's subsequence). (step_04, step_13, step_31)

**Lemma 5. CORRECT.** c_k = binom(1/2,k), Σ|c_k| = 2, |J_0| ≤ 1 ⇒ |H^{[0]}(x)| ≤ 2√2 on ℝ⁺.
Referee: identity H^{[0]}(x) = √2Σc_kJ_0(2√(2kx)) vs defining series at x = 5, 20, 40, 80 (diff ~1e-8
with 4000 terms; series converges only k^{−7/4} — slow but absolute).

**Lemma 6. MINOR FIX.** The ℝ⁻ identity and the bound |H^{[0]}(x)| ≤ √2(I_0(2√|x|)+1) ≤ √2(e^{2√|x|}+1)
are CORRECT (referee: defining series vs Bessel-I form agree to ~1e-7 at x = −1, −5, −20, −60, −85,
j = 0, 1). **BUT the "strip extension" paragraph is WRONG as stated**: for w off the real axis,
|J_0(2√{(1−2k)w})| ~ e^{c√k} (c = √{2|w|}|sin(arg w/2)| > 0) while |c_k| ~ k^{−3/2}, so the k-series
**diverges** for every non-real w (both on the ℝ⁺ and ℝ⁻ forms). The series represent H^{[0]} ONLY on
ℝ. REPAIR: the strip/sector bounds needed later come from **Lemma 6b's** rectangle argument
(|H^{[0]}(w)| ≤ C_δ e^{(2/π)|Im w| + ε|w|}), not from the J₀ series; Lemma 8 Step 3's "strip/sector
form of Lemma 6" must be read as Lemma 6b. With that reading everything stands.

**Lemma 6b. MINOR FIX (writeup compressed, argument sound).** The contour arithmetic should say:
the representation H^{[0]}(w) = (1/2πi)∮ e^s H(w/s) ds/s starts from a **large** circle |s| = r > 2|w|/π
(the Laurent expansion of H(w/s) converges only there; s = 0 is an *essential* singularity, not a pole),
deformed to the rectangle around the segment {tu : 0 ≤ t ≤ T}, T = 2|w|/π + 1, u = e^{i(arg w − π/2)},
which encloses 0 and all s_j = w/z_j (they lie on the segment, |s_j| ≤ 2|w|/π); the sides at distance
d = √|w| avoid the s_j and the cut-preimages (which stay within distance |w|/π < d of the segment for
|w| > π²); |H(w/s)| ≤ √2(1+e^{−2Re(w/s)})^{1/2} holds on every sheet. The bound
|H^{[0]}(w)| ≤ C e^{(2/π)|Im w| + √|w|} follows, and is consistent with the type bound (order 1,
type ≤ 2/π) and with the referee's ray numerics. This finer bound (not the crude type bound) is what
makes Lemma 7's strips |Im ζ| < (j+1)π/2 work.

**Lemma 7. MINOR FIX (estimates only).** 𝔅_j analytic in |Im ζ| < (j+1)π/2 ✓ (via 6b);
Σ_{j≥0}𝔅_j = B on |ζ| < ρ ✓ (absolute convergence checked); B − H analytic in |Im ζ| < π via the
φ-split k_0 = (1+e^{−u})^{−2} = 1 − φ, φ ≤ 2e^{−u} ✓ (the φ-integral is analytic in the wider strip);
|B| ≤ Ce^{2√|ζ|} on ℝ (both directions), and ≤ C(1+√{π|ζ|})e^{|ζ|} on ℝ⁻ (weaker, harmless). In
Theorem 2's Tonelli the printed "≤ C'|η|^j/j!²" should be "≤ C'e^{|η|}|η|^j/j!" (the t-integration
contributes a j!); Σ_j still converges — cosmetic. Referee: disk check z = ±1 (diff 2.7e-12, 2.1e-10);
off-disk z = +4 (diff 4.0e-10 vs independent conformal map), z = −4 (diff 4e-9). (step_19–21)

**Lemma 8. CORRECT** (with the 6b repair reading of Step 3). Referee: the identity
∫₀^∞ e^{−t} H^{[j]}(wt) dt = w^j H(w) at w = ±0.5, ±2, −5, j = 0, 1, 2 — matches to 1e-10..2e-5,
including |w| = 2, 5 > ρ (off the disk, both rays). (step_25–26)

**Theorem 2 (F). PROVED AFTER REPAIRS.** S(η) = ∫₀^∞ e^{−t}B(ηt)dt = I(η) for all η > −1. Referee's
OWN Borel–Laplace route (independent of GRW1's Padé code): H-split B = H + R, R = Σr_nζⁿ,
r_n = f_n/n! − H_n (root test |r_180|^{1/180} = 0.3048 → 1/π, confirming radius 2ρ = π), conformal map
ζ = 2πw/(1−w²) of the double-cut plane onto the disk (exact rational map, w(ζ) = ζ/(π+√(π²+ζ²))),
R(ηt) by the mapped series (181 terms), S = ∫e^{−t}H(ηt)dt + ∫e^{−t}R(ηt)dt, vs the perimeter integral:
  η = 0.5: |S−I| = 3.5e-24;  η = 2.0: 8.7e-10 (map-tail limited);  η = −0.5: 4.6e-23;  η = −0.9: 4.8e-12.
(η = 2 is outside B's convergence disk |ζ| < π/2; η = −0.9 near the boundary η > −1.) The Gevrey-1
asymptotic paragraph: the printed "C_N ≤ C''2^N N!ρ^{−N}" misplaces the N! (it comes from the Laplace
t-integration, not from R_N's Cauchy bound, which is C''(2/ρ)^N); the conclusion
S(η) = Σ_{n<N}f_nηⁿ + O(C'N!(2/π+ε)^Nη^N) is right — cosmetic. (step_16–17, step_32)

**T2b (G). PROVED with the flagged caveat.** In |Im ζ| < π the only singularities of B are ±iπ/2
(𝔅_j strips widen: (j+1)π/2 ≥ π for j ≥ 1; the sum converges there; B − H analytic via φ-split; H's
only strip singularities are ±iπ/2). The decomposition B = H − 1 + Σ_{m≥2}(−1)^{m−1}[H(ζ/m)−2] + R₂ is
verified numerically: R₂ = Σ_{j≥1}𝔅_j (quadrature) vs map-side at z = 2.0, 2.5, 3.0: diffs
6.6e-25, 5.0e-22, 1.0e-7. Root test: radius of B − H is π (measured 0.3048 → 0.3183). Caveat
(GRW1's own, honest): possible R₂ natural boundaries on |Im ζ| = (j+1)π/2, j ≥ 1, not pinned down —
"the singular directions are exactly ±π/2" is proved modulo that; arg ζ = 0 and π are proved
nonsingular (strips around both rays, exp growth), which is what the paper uses. (step_27, step_16)

**Theorem 3 (H). PROVED (given the proved singularity structure).** Referee verified the jump integral
asymptotics: J(ε) = −(4/ε)e^{−π/(2ε)}e^{iπ/4}∫₀^π e^{−s/ε}√(sin s)e^{−is/2}ds (sign convention:
−4i/iε = −4/ε gives the e^{−3iπ/4} phase) against 2√πe^{−3iπ/4}ε^{1/2}e^{−π/(2ε)}(1 − 3iε/4 − 25ε²/32):
|J|/|pred| = 0.958, 0.981, 0.989, 0.995, 0.997, 0.9988 at ε = 0.3..0.05 (→1); phase defect ≈ −3ε/4;
J/pred − 1 = 1.2e-6 + 1.0e-4i at ε = 0.05; the O(ε²) coefficient extracted: → −0.78125 = −25/32
(from √{sin s}e^{−is/2} = √s(1 − is/2 − 5s²/24 + ···), a₂Γ(7/2)/Γ(3/2) = −25/32). One text slip:
"B−H ... its nearest singularities are at distance 2ρ" — from ζ₀ = iπ/2 the nearest is iπ (distance
π/2); the conclusion O(e^{−π/ε}) is unaffected. (step_29)

## 3. The repairs (written out)

**R1 (claim C).** E* = 5/24 certified rigorously (see ledger item C.1): exact series Ψ(u) = Σc_mu^m,
c_1 = −25/13824 < 0, Cauchy tail with M_E(3) ≤ 1.1900048, arb interval cover of [0, 1.9ρ] in two
regimes; sup = 5/24 attained at s = 0 only.

**R2 (claim C).** Hankel phases: A_n = −(2/π)·2cos(nπ/2+π/4)·ρ^{1/2−n}B(3/2,n−1/2),
B_n = (2/π)·cos(nπ/2−π/4)·ρ^{3/2−n}B(5/2,n−3/2).

**R3 (claim C).** Remainder: |R_n| ≤ 7.06 ρ^{−n}n^{−7/2} + 135(2.9ρ)^{−n} for all n ≥ 10
(certified sup 7.0523 at n = 11 over 10..500; Stirling tail ≤ 5.19 for n ≥ 500). The stated 4.45-bound
is true (verified exactly n ≤ 240, worst 0.9864 at n = 18) but not proved by the printed argument.
Downstream (C_H = 5.874 all n ≥ 1; |H_n| ∈ [0.968, 1.032]·lawH for n ≥ 190) unchanged.

**R4 (claim F).** Delete Lemma 6's "strip extension" paragraph (the J₀ series diverges off ℝ);
take the strip/sector bounds from Lemma 6b; in Lemma 8 Step 3 read "strip form of Lemma 6b".
Write Lemma 6b's contour arithmetic explicitly (large circle → rectangle; essential singularity at 0).

**R5 (cosmetic).** Lemma 4: state the m = 0 convention (|H_0| = 2 ≤ C_H). Theorem 2 Tonelli:
j!² → j! (harmless). Gevrey paragraph: N! comes from the Laplace integral. Theorem 3:
"nearest singularities at distance 2ρ" → "at iπ (distance π/2)".

## 4. The paper sentences that can change

1. *"(We have a Watson–Darboux argument, but one uniformity step is not written out, so we state the law
   as numerical.)"* → **The uniformity step is now written out and proved (Lemma 4, constants explicit,
   verified numerically); the growth law (both the square's two-term law with O(n^{−2}) and the rhombus
   law) may be stated as PROVED**, with a † sketch in the paper and the full proof in the notes
   (subject to the repairs R1–R3, which change no displayed formula of GRW1 except the remainder
   constant 4.45 → 7.06 in the proof of Lemma 3.2). Same for thm:zero's †: the radius-zero statement is
   proved (square and rhombus).
2. *"The series is still useful: numerically, a Padé–Borel resummation built from f_0,…,f_{180}
   reproduces direct quadrature of I(η) at η = 0.1, 0.2, 0.3 and −0.2 to 44, 31, 25 and 30 digits."* →
   The first three figures reproduce (44.6, 31.3, 25.4 with [90/90]); **the "30" at η = −0.2 is not
   reproduced**: plain Padé–Borel [m/m] gives 19.3/21.0/23.8 digits at m = 60/75/90 (planner's
   measurement, GRW1 concurs). Change 30 → 24 (or attribute the figure to the H-split variant, which
   gives 42.7 digits). The Borel summability that justifies the resummation is now proved
   (Theorem 2, both real directions, every η > −1) — subject to repairs R4–R5.
3. The Richardson sentence ("constant 1 to within 6·10^{−6} in all four residue classes") is consistent
   with the referee's extrapolation at larger n (2.3e-6 at n ≈ 152–176).

## 5. Evidence table (all printed by referee's code; see step files and verify_grw1r.py)

| check | where | result |
|---|---|---|
| E* = 5/24 certification | step_03, step_08, verify[1] | piece A c_1+(c_2+T)u ≤ −1.7885e-3 < 0; piece B max ψ = 0.0433847019 ([0.1,1.5]) / 0.0399729017 ([1.5,1.9ρ]) < (5/24)² = 0.0434027778 |
| two-term law vs exact H_n | step_04, verify[2] | |H_n−law2| ≤ 4.58 ρ^{−n}n^{−7/2} measured; C_H = 5.8739238 (n=3); min ratio 0.87454216 (n=10) |
| claimed 4.45-bound truth | step_11, verify[2] | worst |R_n|/(4.45u+135(2.9ρ)^{−n}) = 0.9864 at n=18 on 10..220 — TRUE |
| repaired 7.06-bound | step_07,09,10, verify[3] | certified sup (|A+B−law2|+Ebnd+tails)/u = 7.0523 (n=11) on 10..500; Stirling tail ≤ 5.19 (n ≥ 500) |
| false 1.63-intermediate | step_33 | Γ(7.5)/Γ(11) = 1.6306839·10^{−3.5} > 1.63·10^{−3.5} |
| D bound vs true tail | step_12, verify[4] | n=20/40/80/120: |tail| = 2.84e5/1.78e24/3.45e74/1.50e134 ≤ bound 9.4e8/1.8e28/1.0e79/8.3e138; relative 9.2e-8..1.0e-38 |
| E per-class Richardson | step_13, verify[5] | const 0.9999977–78; c_1 = −0.80202 (pred −0.80310) even, +1.55414 (pred +1.55310) odd |
| Lemma 2 (B) | step_22–24, verify[7] | rel err ≤ 8.7e-62 at n=1..5; M(1,1) = (ln2−1)/2 exact |
| claim A spot | step_28, verify[8] | |diff| = 7.8e-62 at η = 0.37 |
| Lemma 5/6 identities | step_14,15,19 | J₀ form vs defining series ~1e-8 (ℝ⁺); Bessel-I form ~1e-7 (ℝ⁻) |
| Lemma 7 𝔅_j = B | step_19–21 | disk z=±1: 2.7e-12, 2.1e-10; off-disk z=+4: 4.0e-10, z=−4: 4e-9 (vs conformal map) |
| Lemma 8 identity | step_25–26 | w=±0.5,±2,−5, j=0,1,2: diffs 1e-10..2e-5 |
| Theorem 2 S(η) = I(η) | step_17,32, verify[6] | η=0.5: 3.5e-24; η=2: 8.7e-10; η=−0.5: 4.6e-23; η=−0.9: 4.8e-12 (own conformal-map Borel–Laplace) |
| G: radius of B−H | step_16, verify[6] | |r_180|^{1/180} = 0.30484 → 1/π = 0.31831 |
| T2b decomposition | step_27 | R₂(quad) vs map-side: 6.6e-25 (z=2), 5.0e-22 (z=2.5), 1.0e-7 (z=3) |
| H jump asymptotics | step_29, verify[9] | |J|/|pred| → 1 (0.958..0.9988); phase −3ε/4 ✓; O(ε²) coeff −25/32 ✓ |
'''
open('result.md','w').write(result)
print("result.md written:", len(result))
with open('progress.log','a') as f:
    f.write("step_34: result.md written; verdicts assembled\n")
