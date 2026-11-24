# Referee report: `thm:pip` (Theorem 14, expansion of π_p) and `prop:ruler` (Proposition 4, two exact rulers)

Paper: `paper/ieee/main.tex`, draft 37, lines 1317–1361 and table rows 2085–2086. Notes read:
`D12_pi_p/result.md`, `RUL1_level2_exact_reductions/result.md` (+ `verify_rul1.py` read, not run or imported),
`D12c_ruler_family/result.md`. All code here is my own (`r1`–`r6`); logs `r*_log.txt`.

## Verdicts

| Result | Verdict |
|---|---|
| `thm:pip` as stated (N_0..N_3, O(p^-5)) | **PROVED** |
| Text claim "N_4, N_5 are proved too" (and table row "N_0,…,N_5") | **PROVED AFTER REPAIRS**: one missing paragraph (all-orders remainder); the Euler-sum reductions are exact, and I re-proved them with my own code |
| `prop:ruler` | **PROVED** (the proof sentence needs a small fix) |

No printed constant rests on an integer-relation fit. D12 found N_2..N_5 and the J values by PSLQ first. The proof of N_0..N_3 is by hand. N_4 and N_5 are exact: RUL1's certificates, and independently my own DS+duality elimination (r2).

## (a) Statement checks

| Test | Script | Result |
|---|---|---|
| π_p from geometry (x-parametrisation, ℓ^p length) vs the split form 4·2^{-1/p}+E(p), p = 3, 5, 10, 100, 10^3, 10^4 | r3 | agree to 80 digits (working precision) |
| π_2 = π; duality π_p = π_{p/(p-1)}, p = 3, 5, 10; minimum at 2 | r3 | 80 digits; minimum confirmed |
| N_0..N_5 closed forms vs the term-by-term coefficients (J integrals by quadrature) | r1 | 165–170 digits |
| E(p) coefficients π²/6 and π²/6+ζ(3)−π²log2/6 | r1 | 170 digits |
| (π_p/4 − series through N_3)·p^5 → N_4 (p = 10^2..10^5) | r1 | 0.34863, 0.349584, 0.3496707, 0.3496793 → N_4 = 0.3496802; the remainder is O(p^-5) |
| Direct π_p at p = 10^15 minus series through N_8 | r1 | 4.70029e-150 = N_9 p^-10. N_4 confirmed to 74 digits and N_5 to 59 digits **from direct values** |
| Remainder vs next term, p = 10^3..10^15, K = 3, 5 | r1 | residual = next term to 4–6 digits at every p |
| D12 table N_6, N_7, N_8 | r1 | 40 digits (the table's precision) |
| Euler sum Σ(−1)^{n+1}H_n/n² = 5ζ(3)/8; ∫log²(1+v)/v = ζ(3)/4; J_{1,k}=k!η(k+2) | r6, r1 | 60, 60 and 160+ digits |
| L_1(p) = 8 (p = 1.5, 2, 3, 7, 50) | r3 | 80 digits |
| L_∞(p) = 8·2^{-1/p} | r3 | exact by hand; q = 10^6 proxy differs by 10^-12 = O(1/q), as expected |
| L_q geometric (x-form) vs the paper's t-integral, 5 pairs (p, q) | r3 | 28–57 digits (raw t-quadrature is the limit) |
| Fixed q = 1.5, 3, 4: 200 Taylor coefficients of F_q from a Cauchy integral on \|ε\| = 0.85; partial sums vs L_q(p) at p = 3, 1.5, 1.25 | r4 | 28–34 digits at 34-digit working precision; errors fall steadily with K at p = 1.25 |
| Growth of N_k, k ≤ 30 (my own numerical series engine, agrees with r1 to 68 digits for k ≤ 8) | r5 | see (c)4: true for odd k only |

## (b) Proof checks

**thm:pip, N_0..N_3 (D12 §2).** I checked every step by hand.
- Lemma 1 (split): exact. 4h∫_0^1 r^{h-1}(1+r)^{-1-h}dr = 4·2^{-h}, since the integrand is d/dr (r/(1+r))^h.
- Lemma 2 (substitution): exact. I re-derived it: 1/(p−1) = h/(1−h) and r^{p−1} = e^{−s}.
- Lemma 3(iii): correct. log w = −log2 + h(s/2 − log2) + O((s+s²)h²) for s ≤ h^{-1/2}, and hs/2 ≤ √h/2 keeps the exponentiation uniform.
- Proposition: correct. The tails beyond h^{-1/2} are O(e^{-h^{-1/2}}) because L ≤ e^{-s}. The algebra gives R = 2J_{1,0}h³ + (J_{1,1} − 2log2·J_{1,0} + 2J_{1,0} + J_{2,0})h⁴ + O(h⁵); r1's sympy expansion gives the same coefficients.
- The J evaluations are proved by hand: the J_{1,k} series, and J_{2,0} via the classical Euler sum (H_{n−1} = H_n − 1/n gives −(5/4)ζ(3) + (3/2)ζ(3) = ζ(3)/4).
- There is no differentiation under the integral sign. The expansion is controlled by explicit bounds.

**N_4, N_5.**
1. *All-orders remainder: not written.* D12 says "same machinery" and labels the formula N_k = (−log2)^{k+1}/(k+1)! + Σ∫[h^m]Ψ as PROVED, but §2 only does orders h³ and h⁴. The gap is routine and closes as follows. Write w = W(h, sh) with W(x, y) = e^{−xy/(1−x)}(1+e^{−y/(1−x)})^{−1−x}. W is analytic for |x| < 1/2 and |y| < π/2. Cauchy estimates bound its Taylor remainder of total degree > m by C·(h(1+s))^{m+1} for s ≤ h^{-1/2}. That is O(h^{m+1}) after integration against L. The tails and the expansion of e^{hL} − 1 are handled as at order 4.
2. *Euler-sum reductions: exact.* My r1 derivation shows which J's are needed:
   - N_4 needs J_{2,1}/8 + J_{3,0}/12;
   - N_5 needs J_{3,1}/24 + J_{4,0}/48, plus J_{2,1}, J_{3,0}, J_{2,0} and J_{1,k}.

   **J_{2,2} is not used.** RUL1 lists it among "the five", and D12 lists it too. The five nontrivial sums actually used are J_{2,0}, J_{2,1}, J_{3,0}, J_{3,1}, J_{4,0}.
3. *Independent exact re-proof (r2).* I wrote my own engine for level-2 iterated integrals G(w;1) with two relation families:
   - finite double shuffle on convergent words;
   - duality under t ↦ (1−t)/(1+t). I re-derived the reversal sign (−1)^n and τ(0) = (1)−(−1), τ(1) = (0)−(−1), τ(−1) = −(−1).

   Results:
   - The ranks give codimensions 1, 2, 3, 5, 8 at weights 1–5, the Deligne–Goncharov bound.
   - Control: ζ(2̄,1) = ζ(3)/8 is proved.
   - For each of k = 0..5, N_k (as an exact rational combination of log2 and the J words) minus the printed closed form lies **exactly** in the Q-span of the relations, weight by weight. The weight-6 part (log⁶2) cancels identically.
   - Four perturbed closed forms (π⁴, ζ(5), ζ(3)π², ζ(3) changed by small rationals) are correctly *not* in the span.

   This route never uses Li_4(1/2) or Li_5(1/2), so RUL1's "70-digit checked" word expansions of Li_n(1/2) are not needed for N_4, N_5. The individual J closed forms (which do contain Li_n(1/2)) I checked only to 160+ digits. They are not needed for the theorem.
4. RUL1's certificates are re-verifiable in principle: they are exact rational identities over explicit relation instances. I did not run `verify_rul1.py`. I reproduced the conclusion with independent code instead.

**prop:ruler.**
- L_1 = 8 and L_∞ = 8·2^{-1/p} for p ≥ 1: correct by hand. On t ≤ 1/2, |dx/dt| ≥ |dy/dt| because the exponent 1/p − 1 ≤ 0. At p = 1 both formulas give the right values (8 and 4).
- Fixed q ∈ [1, ∞): put z = (t/(1−t))^{q(1−ε)}, so that |z| ≤ 1 for Re ε < 1. The integrand of F_q = ∫(a+b−(a^q+b^q)^{1/q}) is b − a·((1+z)^{1/q} − 1) = b − (az)·g(z), with g(z) = ((1+z)^{1/q} − 1)/z bounded on the closed unit disc. Also |az| = |b|^q|a|^{1−q} ≤ |b| ≤ 2^{1−Re ε}, **using q ≥ 1**. So F_q is holomorphic on Re ε < 1, the series converges for |ε| < 1, and its sum is L_q. Correct.
- The paper's one-line proof ("works with (1+z)^{1/q} in place of √(1+z)") is too quick. The proof of `lem:holo` uses the identity b − b²/(a(1+√(1+z))) and z = (b/a)², which do not carry over verbatim. The q-version needs z = (b/a)^q and the bound above. See the repair below.

## (c) Sentences that claim more than is established (line numbers of main.tex)

1. **L1351** "The next two coefficients … they are proved too: the five alternating Euler sums they need reduce exactly …". This is true once the all-orders remainder is written. Suggested: "The next two coefficients bring in π⁴ and ζ(5). They are proved too (in our notes): the same bounds, carried to order p^{-7}, reduce them to five alternating Euler sums, and double shuffle and duality relations reduce those exactly to Li_5(½), Li_4(½), ζ(5), ζ(3), π² and log 2. The polylogarithms cancel." Also **print N_4 and N_5** (formulas below), or change table row L2086 to "N_0,…,N_3 (N_4, N_5 in notes)". At present the row claims N_0..N_5, but the theorem states only N_0..N_3.
2. **L1351** "it enters π_p at order p^{-4}, and N_3 contains ζ(3)/4 in the basis {ζ(3), π²log2, log³2}". N_3 is not of pure weight, and log³2 does not occur in it. "Enters at p^{-4}" presupposes that ζ(3) is not a Q-combination of π² and log³2 (unknown). Suggested: "in our closed forms ζ(3) first appears in N_3, whose weight-three part is ζ(3)/4 − π²log2/24."
3. **L1351** "Its constants are now log2, π², ζ(3), ζ(5): there is no √2 and no log(1+√2)." Suggested: "N_0,…,N_5 are polynomials with rational coefficients in log 2, π², ζ(3) and ζ(5); no √2 or log(1+√2) is needed."
4. **L1351** "this series diverges: numerically |N_k| grows like 2(k−2)!/π^{k−1}". Divergence is not proved (D12 labels it NUMERICAL). The growth law holds **only for odd k**. My r5 gives N_kπ^{k−1}/(2(k−2)!) as follows:
   - odd k: +1.283 (k=21), +1.219 (k=25), +1.179 (k=29), with signs following the period-four pattern;
   - even k: +0.509 (k=22), +0.430 (k=26), +0.375 (k=30). Here the normalised value is about 11/k, so |N_k| is smaller by a factor of about 11/k.

   Suggested: "and unlike P(p), this series appears to diverge: numerically (k ≤ 40), |N_k| ≈ 2(k−2)!/π^{k−1} for odd k, the even-k coefficients are smaller by a factor of order 1/k, and the signs have period four."
5. **L1353, L1361** "Why does π_p diverge when P(p) converges?" / "this is why π_p diverges". π_p is a finite number. Suggested: "Why does the 1/p series of π_p diverge when that of P(p) converges?" and "this is why the series of π_p diverges while that of P(p) converges; we have not proved this."
6. **L1359** proof sentence for prop:ruler. Suggested: "…and for fixed q ≥ 1 the proof of Lemma~\ref{lem:holo} works with z = (b_ε/a_ε)^q = (t/(1−t))^{q(1−ε)} and (1+z)^{1/q} in place of √(1+z): the integrand is b_ε − a_εz·g(z), with g(z) = ((1+z)^{1/q}−1)/z bounded on |z| ≤ 1 and |a_εz| = |b_ε|^q|a_ε|^{1−q} ≤ |b_ε| for q ≥ 1."
7. **L1361** "So neither extreme ruler sees ζ(3)." This is acceptable as a statement about formulas. Clearer: "The coefficients of both extreme rulers are rational multiples of powers of log 2."
8. **L1325** "with an explicit integral E(p)": the integral is never displayed. Suggest displaying it (verified here to 80 digits):
   E(p) = (4/(p(p−1))) ∫_0^∞ e^{−s/(p(p−1))} (1+e^{−s/(p−1)})^{−1−1/p} [(1+e^{−s})^{1/p} − 1] ds.

Checked and fine: the definitions of π_p and L_q (re-derived from ℓ^q arclength), L_2 = P, L_p = 2π_p, the literature facts (π_1 = π_∞ = 4, π_2 = π, minimum at 2, duality), the integral and Euler-sum formulas in the proof sketch, and "for every p ≥ 1" in prop:ruler (exact for all p ≥ 1).

**N_4, N_5 (verified exactly in r2 and to 165+ digits in r1):**
- N_4 = −log⁵2/120 + (π²/24)(1 − log2 + log²2/2) + (ζ(3)/4)(1 − log2) + π⁴/2880
- N_5 = log⁶2/720 + (π²/24)(1 − log2 + log²2/2 − log³2/6) + (ζ(3)/4)(1 − log2 + log²2/2) + (π⁴/2880)(1 − log2) − 7π⁴/1920 − ζ(3)π²/48 + ζ(5)/8

## Files
`r1_pip_series.py` (direct values, sympy term-by-term expansion, J quadrature, remainder checks), `r2_exact_words.py` (exact DS+duality proof of N_0..N_5), `r3_direct.py` (geometric π_p, L_1, L_∞, L_q), `r4_ruler.py` (fixed-q convergence), `r5_growth.py` (N_k to k = 30), `r6_euler.py` (Euler sum). Logs: `r*_log.txt`. Data: `r1_Nsym.pkl`, `r1_Nnum.pkl`.
