# M5AR: referee report on M5A ("M₅ at every angle")

An AI re-check with its own code. It is not peer review. The claim is in `new_results/M5A_M5_every_angle/`. I changed no existing file. Everything below is in `new_results/M5AR_referee_M5_every_angle/`. The claim's files were copied into `claim_copy/` and run only from there. The rerun output is `claim_copy/s03_rerun.out`.

## Verdict: CONFIRMED

The mathematics is confirmed, with no mathematical repairs needed. Several sentences of `result.md` claim more than was shown, or give an incomplete reason (findings 1–4). They need rewording, not new mathematics.

- **The object is right.** I computed M₅(α) from the paper's definition with my own code. It does not use μ₅,ₙ and has no Cauchy integral. It agrees with the claimed closed form, evaluated by my own HPL evaluator, to within 1.3·10⁻⁹⁸ (relative ≤ 5.3·10⁻⁹¹) at eight angles. The angles run from α = 0.2 (c = 0.995) to α = 3.1 (c = 0.0208). At the square corner the definition agrees with the stored 400-digit M₅ to 5.4·10⁻⁷⁸ at 80 dps.
- **"Proved for all n" is sound.** The proof is a word-by-word identity in Q(n) of Z-sum forms. Word-by-word equality is a sufficient condition, so it is valid for every n ≥ 1 whatever the words are. The rerun is identical. All 12 planted errors were caught, and the unmutated control passed.
- **Independent exact check.** My own derivation of μ₅,ₙ from the definition (Beta integrals and polygammas) matches every part, as exact rationals, at every Laurent coefficient c^N with −6 ≤ N ≤ 301. That is all n ≤ 150.

## Findings (most serious first)

1. **The reason given for "the π⁴ bracket is no longer elementary" is incomplete.** The statement is true. (`result.md` l.114)
   - The text argues only that D appears with a nonzero coefficient. That does not exclude D being elementary modulo the other functions.
   - **Fix:** add one sentence. By the linear independence of harmonic polylogarithms over C(c) (Deneufchâtel–Duchamp–Hoang Ngoc Minh–Solomon, 2011), an identity bracket = r₁ + r₂·artanh c with rational r₁, r₂ would be a nontrivial C(c)-linear relation among H₀,₁, H₀,₋₁, H₋₁,₁, H₁,₋₁ (which occur only in D), H_{±1} and 1. That is impossible because the coefficient (1−c²)(3c²−7)/(1920c³) is not zero.
   - The same argument shows the bracket is not a polynomial in logarithms either. At the symbol level, the antisymmetric part of S(D) is −(1−c)∧c + (1+c)∧c + 2(1−c)∧(1+c) ≠ 0.
   - I verified the explicit π⁴ bracket exactly through c^301 (`r03` D).

2. **"This is what forces e ∈ {−4,…,1}" (l.49) is not shown.** The denominators 2n+3 and n+2 of μ₅,ₙ make the range plausible, but the e-range was an assumption of the fit.
   - **Fix:** either delete "forces", or replace it with the correct argument. By the same linear-independence theorem, the rational function multiplying each H_w(c) is unique. So the c⁻⁴ prefactors in the proved representation are intrinsic, not an artefact of the basis, for representations by H_w(±c) with words in 0, ±1 not ending in 0.
   - The same remark upgrades the "unique (rank r/r)" statements in §3 and l.67 from "unique in this finite basis" to "unique as rational-function coefficients of HPLs".

3. **"Field rule at k = 5 (PROVED) … as for k ≤ 4" (l.112) understates what is known and overstates what is new.** The paper's proof sketch of `conj:field` (main.tex l.1488) already says: "For k≤7 the denominators are powers of c alone".
   - **Fix:** say that the proved k = 5 representation has prefactors c^e, −4 ≤ e ≤ 1, consistent with that statement. By finding 2 these prefactors are canonical. "The range grows by one" is then a correct statement about these canonical representations.

4. **The ζ(4)/ζ(2)² convention is stated (l.15) but not defined.**
   - **Fix:** say "formal coefficient in the sense of Definition g6:def:formal: each ψ^(m−1) at an integer contributes ζ(m) as a symbol, and products come from products". This is the convention of Theorem `thm:zladder`'s kernel log Γ(1+ξ) + log Γ(1+η) − log Γ(1+u) = Σ (−1)^m ζ(m)/m·(ξ^m + η^m − u^m).
   - Also say which statements depend on it. The convention matters for the rung [ζ4]M₅ = Ψ + G, the values of G and Ψ+G separately, and "[ζ(4)]M₅ contains D". It does not matter for the π⁴ bracket π⁴(Ψ/90 + 7G/180), nor for the ζ(5) and ζ(2)ζ(3) parts.
   - My independent derivation computes the split from Beta functions at integer points. It gives exactly the same split as the claim for every n ≤ 150. So the split is robust between these two natural routes, but it remains a convention.

5. **The top-part paragraph (l.115–122) mixes statuses.** The three formulas for [ζ2^e]P₂ₑ,ₙ come from the claim's symbolic μ-side (exact for all n, by computer algebra, if that code is right at k = 6). I confirmed [ζ2³]μ₆,ₙ = n(5n−1)/(3(2n−1)(2n+1)(2n+3)(2n+5)) independently for n ≤ 80.
   - The sentence "By the ladder, every top part is a rational combination of E_e" is asserted without proof.
   - **Fix:** label it as a remark, or give the one-line weight count from the kernel.

6. **Minor.**
   - l.129 "No referee has checked this work" should be updated.
   - In s01, the μ-side is tested against `g6_common` only for n ≤ 8. That test is not independent of the derivation's lemmas. This report now supplies n ≤ 150 from a different derivation.
   - The proof's μ-side at k = 5 still rests on the general-k lemmas of `sec_g6_flat` (Lemmas g6:lem:En and g6:lem:A), which I did not re-read (see "Not checked").

No error was found in any formula, coefficient, count or number in `result.md`. These all match my output: the term counts (310 = 1+8+27+84+190; 60 = 3+6+19+32; 13; 22; 5; 5), the denominators of R₅ (1, 2, 4), G and the [ζ5] and [ζ2ζ3] closed forms, the Ψ closed form, the s03 formulas for [ζ2²]μ₅,ₙ and [ζ4]μ₅,ₙ (with Z₋₁(2n)), the "−3/c³+2/c+c" remark, and the six M₅ values in §5.

## What I verified and how

The scripts are in this folder and run with `../../../.venv/bin/python`. None of r01–r04 or r06 imports `m5fit`, `g8_allN` or `g6_common`.

**r01_mu_exact.py: the μ_{k,n} side (item 3).** Exact, k ≤ 5, n ≤ 150, about 1 min.
- *Starting point.* The length of x e₁ + y e₂ with x = t^u, y = (1−t)^u is 2 + 2uM. Using (a+b)² − 4c²ab, this gives μ_n(u) = ½·binom(½,n)·(−4)^n·J_n(u), with J_n = ∫₀^∞ g^{(n−1)(1−u)}(1+g)^{−1−u}(1+g^{1−u})^{1−2n} dg.
- *Expansion.* I expanded (1+g^{1−u})^{1−2n} in (g^{−u}−1)^i with i ≤ k. Each piece is then a Beta integral Γ(n+i+u(1−n−j))Γ(n+u(n+j))/Γ(2n+i+u).
- *Zeta symbols.* I used ψ^(m−1)(N) = (−1)^m(m−1)!(ζ(m) − H^(m)_{N−1}). γ cancels exactly.
- *Sanity checks.* μ₀,ₙ = −1/(n(2n−1)), [ζ2]μ₂,ₙ = 1/(4n²−1) and [ζ3]μ₃,ₙ = −1/(4n²−1) all hold (`thm:flat`). The value of M₀ also checks numerically.

**r03_series_compare.py: exact coefficient comparison (items 3 and 4).**
- *Method.* My own Taylor coefficients of H_w come from the definition. I used my own reading of the pair, c^e H_w(c) + (−c)^e H_w(−c), not the claim's sign s. I compared all coefficients c^N, −6 ≤ N ≤ 301, exactly.
- *Results.* All 26 checks are True:
  - each of the 7 fitted files against my μ₅;
  - exactly 7 zeta monomials occur in μ₅,ₙ;
  - [ζ5]μ₅,ₙ = −3/((2n−1)(2n+1)(2n+3)) and [ζ2ζ3]μ₅,ₙ = 2n/((2n−1)(2n+1)(2n+3)), for n ≤ 150;
  - [ζ2²]μ₄,ₙ = −n/((2n−1)(2n+1)(2n+3));
  - F₂ (closed form) = [ζ2²]M₄, Φ = [ζ2]M₂, and G8's Ψ and Ψ₂ files = [ζ2]M₃ and [ζ2]M₄;
  - the theorem as stated, with these: −Ψ₂ = [ζ3]M₅, Ψ + G = [ζ4]M₅, −Φ − 2F₂ = [ζ5]M₅, −2F₂ = [ζ2ζ3]M₅;
  - the paper's example rung [ζ5]M₅ = −[ζ2]M₂ − 2[ζ2²]M₄, and [ζ3]M_{k+1} = −[ζ2]M_k for k = 2, 3, 4, on my μ;
  - the closed forms of G, [ζ5], [ζ2ζ3], Ψ and the π⁴ bracket in result.md;
  - the paper's `thm:M4alpha` π⁴ bracket.
- *Bug in my first run.* It reported mismatches, but the cause was in my code: `(-1)**e` is a float in Python for e < 0. I fixed it (integer parity), and the rerun is the output above.

**r02_definition.py: M_k(α) from the paper's integral (item 1).**
- *Method.* At each node g I take the exact truncated u-series of the integrand. I use log1p/expm1 so nothing cancels near g = 0. I integrate the u^k coefficient with tanh-sinh on [0, ⅛, ½, ⅞, 1].
- *Checks.* M₀(π/2) is exact to 30 digits. M₅(π/2) agrees with `../M5a_next_coefficient/M5_400.txt` to 1.7·10⁻⁵⁷ at 60 dps and to 5.4·10⁻⁷⁸ at 80 dps.

**r04_numeric_compare.py: closed form against the definition (item 1).**
- *HPL evaluator.* My own: fixed-point 400-bit integer Taylor sums of H_w(±c), depth-first over suffixes, from the recursion of the definition.
- *Combination.* The seven files, with real ζ values, against r02 at 80 and 100 dps.
- *Angles.* α = 0.2, 0.5, π/3, π/2, 2π/3, 2.5, 3.0, 3.1.
- *Results.* At every angle |closed − def₁₀₀| ≤ 2.9·10⁻⁹⁷ and |def₈₀ − def₁₀₀| = 5.4·10⁻⁷⁸. This is a constant quadrature error, the same at every angle. Relative differences are at most 5.3·10⁻⁹¹ (at α = 3.1, where M₅ = −2.48·10⁻⁸).

**Item 2: rerunning and testing the all-n proof.**
- *Rerun.* `claim_copy/s03_rerun.out` is identical to `s03_prove.out` apart from timings. All 7 parts pass the even, odd and low-order checks, all 4 ladder rungs hold on the μ side, and all 4 fit-vs-ladder comparisons hold.
- *What the finite computation is.* For each zeta monomial z there are two finite dicts {Z-word: element of Q(n)}:
  - the HPL side at N = 2n, built from h_w(2n+δ) with δ = −e ∈ {−1,…,4};
  - the μ side, P₅,ₙ from Definition g6:def:formal for symbolic n.

  Their difference is the empty dict, and the dict at N = 2n+1 is empty. The coefficients of c^N for −6 ≤ N ≤ 1 are 0 in exact rationals.
- *Why that covers every n.* Equality in Q(n) of every word coefficient implies numeric equality at every n ≥ 1 where no denominator vanishes. "Bounded depth/indices, finite set of rational functions" is automatic for a finite dict. I checked the hypotheses for k = 5, not just inherited them:
  - **Denominators.** On the μ side they are 2n−1, n, 2n+1, n+1, 2n+3 and n+2 (s01). On the HPL side they are 2n+t with −1 ≤ t ≤ 5. None vanishes for n ≥ 1, and all poles are at N ≤ 1, as Lemma g8:prop:crit requires.
  - **Index ranges.** h_w(N) is used only at N = 2n+δ ≥ 2n−1 ≥ 1. The negative shifts use t ∈ {−1, 0}, so M = 2n+t ≥ 1.
  - **Low orders.** e ≥ −4 and H_w = O(c), so the lowest power is c⁻³. Constants sit at N = 0 and −2. All of these lie inside the low-order window −6…1.
  - **Odd coefficients.** N = 2n+1 ≥ 3 is covered by the odd identity, and N = 1, −1, −3 by the low-order window.
- *Planted errors.* `r05_mutants.py` (output `r05_mutants.out`) runs the claim's own checker on modified copies. The unmutated R₅ passes, and **12 of 12 mutants are caught**:
  1. R₅, coefficient of `c^-3 He[0,1,0,1,0,1]` +1/4: even and low;
  2. Ψ₃, coefficient of `c^-3 He[0,1,0,1]` doubled: even and low;
  3. Ψ₃, constant ½ → ⅓: low only;
  4. G, `c^-3 He[0,1]` −3/32 → −2/32: even and low;
  5. R₅, dropped `c^-2 He[1,0,1]`: even and low;
  6. Ψ₃, `c^0 He[1]` → `c^-1 He[1]`: even and low;
  7. R₅, `c^1 He[0,0,0,0,0,1]` → `c^-1` (same parity): even and low;
  8. R₅, argument sign swapped, `c^-4 He[1,0,1,0,1]` → `He[-1,0,-1,0,-1]`: even and low;
  9. G, `c^-4 He[1]` → `He[-1]`: even and low;
  10. [ζ3], an even pair replaced by its plain half: even, odd and low;
  11. [ζ5], coefficient +1/16: even and low;
  12. G, constant c⁻² −¾ → −½: low only.

  As a further check, comparing R₅ against the wrong target (the [ζ2] part) fails, as it should. A sign swap of the argument changes a pair only when (−1)^{e+ν(w)} = −1, so I chose such terms.

**r06_side_checks.py: side claims (item 5).**
- **(a)** The s03 formulas for [ζ2²]μ₅,ₙ and [ζ4]μ₅,ₙ, with Z₋₁(2n), equal my μ for n ≤ 150.
- **(b)** The top part [ζ2³]μ₆,ₙ holds for n ≤ 80 with my μ at k = 6.
- **(c)** Ranks: I computed the rank modulo the prime p = 2⁶¹−1 = 2305843009213693951 of the coefficient map of each fit basis, over all orders c^N (odd orders included), with my own HPL coefficients:

  | Basis | Rank |
  |---|---|
  | W = 1 | 8 / 8 |
  | W = 2 | 26 / 26 |
  | W = 3 | 80 / 80 |
  | W = 4 | 242 / 242 |
  | W = 6 (R₅) | 2186 / 2186 |

  Full rank mod p implies full rank over Q, so the "unique in the basis" statements are proved. They also follow from HPL independence (finding 2).

**Side claims, my status.**
- *Denominators are powers of c at k = 5.* Proved, as a property of the proved representation, which is canonical by finding 2. It is not new: the paper already says this for k ≤ 7.
- *π⁴ bracket not elementary.* True, with the one-line independence argument of finding 1.
- *Uniqueness.* Proved (ranks, reproduced).
- *[ζ2³]P₆,ₙ formula.* Exact for all n by the claim's symbolic computation. I checked it independently for n ≤ 80.

## Proposed statement for the paper (in the style of `thm:M4alpha`)

```latex
\begin{theorem}[$M_5$ at every angle$^\dagger$]\label{thm:M5alpha}
For $0<\alpha<\pi$, with $c=\cos\frac\alpha2$, $\Phi$ as in Theorem~\ref{thm:flat},
$\Psi$ and $\Psi_2$ the $\zeta(2)$ parts of $M_3$ and $M_4$ (Theorem~\ref{thm:M4alpha}), and
$F_2=\bigl(-\frac c{16}-\frac1{8c}+\frac3{16c^3}\bigr)\operatorname{artanh}c+\frac1{16}-\frac3{16c^2}$
the $\zeta(2)^2$ part of $M_4$,
\[
\begin{split}
M_5(\alpha)={}&R_5(c)+\zeta(2)\Psi_3(c)-\zeta(3)\Psi_2(c)+\zeta(4)\bigl[\Psi(c)+G(c)\bigr]+\zeta(2)^2G(c)\\
&-\zeta(5)\bigl[\Phi(c)+2F_2(c)\bigr]-2\zeta(2)\zeta(3)F_2(c),
\end{split}
\]
where, with $L=\log(1-c^2)$, $A=\operatorname{artanh}c$ and
$D=\Li_2(c)-\Li_2(-c)-H_{-1,1}(c)+H_{1,-1}(c)$,
\[
G=\frac7{16}-\frac3{4c^2}-\frac{16-19c^2+c^4}{32c^4}\,L+\frac{7-5c^2}{16c^3}\,A
+\frac{(c^2+3)(c^2-1)}{32c^3}\,D,
\]
$R_5$ is a rational combination of $310$ terms $c^eH_w(\pm c)$ with harmonic polylogarithms of
weights two to six and $-4\le e\le1$, and $\Psi_3$ is a rational combination of $60$ such terms of
weights one to four and the constant $\frac12$ (in our notes). The zeta values are formal symbols as
in Theorem~\ref{thm:zladder} ($\zeta(4)$ and $\zeta(2)^2$ are different monomials); the formula with
real zeta values follows.
\end{theorem}

\begin{proof}[Proof idea]
As for Theorem~\ref{thm:M4alpha}: for each zeta monomial, each Taylor coefficient of the two sides is
an identity of alternating harmonic sums in $n$, which computer algebra proves for every $n$ (the odd
and the finitely many low coefficients vanish). The $\zeta(3)$, $\zeta(4)$, $\zeta(5)$ and
$\zeta(2)\zeta(3)$ parts are those of Theorem~\ref{thm:zladder}. As checks, an independent
computation from the definition agrees exactly through $c^{300}$, and to $10^{-90}$ at eight angles
from $\alpha=0.2$ to $\alpha=3.1$.
\end{proof}
```

Optional sentence after it, which is exactly true: "With real zeta values the $\pi^4$ part is $\pi^4(\Psi/90+7G/180)$. Unlike at $k=4$ it is not elementary: it contains $D$ with coefficient $(1-c^2)(3c^2-7)/(1920c^3)\ne0$, and harmonic polylogarithms are linearly independent over $\mathbb Q(c)$."

Status words I would accept:
- for the theorem: "Proved$^\dagger$ (computer-assisted; refereed once by an independent AI re-check)";
- for the term lists: "data, unique as coefficients of harmonic polylogarithms";
- not "PROVED" in capitals for the top-part paragraph's "by the ladder…" sentence.

## What I did not check

- I did not re-read the written proofs of the general-k lemmas behind `g8_allN.mu_side` (sec_g6_flat: Lemmas g6:lem:En and g6:lem:A, Definition g6:def:formal). That the k = 5 μ-side is right for **every** n rests on them and on that code. My independent agreement covers n ≤ 150 for all seven monomials, and n ≤ 80 at k = 6 for the top part.
- I did not audit `g8_zsum.py`'s shift and stuffle rules line by line. Their k = 5 use is covered by my exact coefficient comparison for n ≤ 150, and by the claim's own 9912-check selftest (s04, not rerun by me).
- Not rerun: s00 (fit self-test), s02 (fits), s04, s05, s06, s07. Their substance is re-derived above with my own code. The exception is the k = 6 symbolic computation, which I checked only for n ≤ 80.
- I did not attempt depth, symbols or cobrackets of R₅. Those are listed as not done in result.md.
- My Theorem `thm:zladder` check is limited: the four k = 5 rungs (exactly, n ≤ 150, and analytically from the kernel's log Γ expansion) and [ζ3]M_{k+1} = −[ζ2]M_k for k ≤ 4. I did not check the general theorem.
