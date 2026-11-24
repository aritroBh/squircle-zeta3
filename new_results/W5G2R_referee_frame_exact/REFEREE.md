# W5G2R: referee report on W5G2 (exact certificates for the 74 conjugated frame relations)

Scope: `../W5G2_frame_exact/result.md`, its scripts t01–t06, `g2common.py`, `g2lin.py`, `dixon2.py`, `w5c2_*.py`,
its logs, `cert_summary.pkl`, `w4cert_summary.pkl`; background W5G1, W5G1R, W5C2, paper/ieee/main.tex (read only).
I did not rerun the 5-hour pipeline. Everything below marked "own" is my code in this folder:
`rc.py` (own word enumeration, Galois action, shuffle, numerical evaluator), `relim.py` (own dense elimination mod p),
`r01`–`r04` (logs `r0N.log`). Python: project `.venv`.

## Verdicts

| claim | verdict |
|---|---|
| G2: each of the 74 vectors u.f_i equals (proved frame rows) + (literal std/lift rows and u=7 conjugates of lift rows) + e ш W, with W in the proved weight-4 span | **CONFIRMED** as an existence statement (independently, mod two primes, all 74, both conjugation sectors). The *exact* rational identities are confirmed by code audit and W5G2's logs only; their large intermediates were deleted and I did not recompute them exactly. |
| (H): every u.f_i is a relation | **CONFIRMED**, under the standing hypotheses listed below (same as Theorem W4R) |
| Corollary W_5 = N_5 + R_5 | **CONFIRMED, conditional** on (c) W5C1's vector D and (d) W5G1's exact split D = J + e ш Z + X, besides the standing hypotheses |
| "No Galois argument is used; only u = 7" | **CONFIRMED** (code audit of t04 `literalize`/`verify`) |
| frame rows proved in domain (t06) | **CONFIRMED** (own domain check, r04), with W4R's caveat on the conversion code |
| result.md's paper sentence | **CONFIRMED WITH REPAIRS** (it omits that W_5 also needs the exact conversion of W_5, N_5, R_5 and W5G1's split; replacement below) |

## 1. The logic chain

From "the printed identities hold exactly" to "each v = u.f_i is a relation among the numbers H(w)".
Conventions (checked against my evaluator): letter −1 = dt/t, letter j = z^j dt/(1−z^j t), z = e^{2πi/8};
H(a_1…a_n) = ∫_{1>t_1>…>t_n>0}; admissible = a_1 ≠ 0, a_n ≠ −1; u acts by j ↦ uj mod 8 fixing 0, −1.

| # | link | status |
|---|---|---|
| L1 | t04 identity v = Σ a_j F_j + Σ c_r r + Σ c'_r (7.r) + e ш W holds on all 46656 words | Code audited: `verify` builds the right side from rows loaded independently (`validated_w5.pkl`, `rels_frame234.pkl`, W5C1 `rels_*.pkl`, `E1`), and compares with `target_vec` = u-image of the stored frame row, on the union of supports, after one common denominator. The solver supplies only numbers (a_j, c_r, W). Not a tautology *provided* W is constrained, which is L6. Exact result taken from t04.log (0 mismatches, 74/74). Independently: existence confirmed mod 2 primes (§3). |
| L2 | `r` ranges only over literal rows; u enters only as u = 7 on lift rows | Confirmed by reading `literalize`: std images are mapped back to literal std rows (`STDK` lookup, asserted); lift images under 3, 5 are replaced by literal lifts of the decomposed base row plus an e-part that goes into W; I-rows go into W; only keys (fam, idx, g) with g ∈ {1, 7} reach `verify`, and g = 7 only for lift families. No literal row is a u = 3 or 5 image of a lift row. Frame basis labels all have u = 1 (`cert_summary`: 4402 framew5 + 966 frame234 terms, all u = 1). |
| L3 | std rows are relations | `make_relations(5, {})`: regularised extended double shuffle reg(u*v − u ш v), u ∈ A^1 (may start with 0), v admissible, and Zhao's regularised distribution, d = 2, 4, 8; final rows on admissible words. Theorems (Racinet; IKZ-type extended double shuffle at level N; Zhao), *plus correctness of relations.py's implementation*, exactly as in Theorem W4R. Own numerics: all 94762 rows vanish, max relative residual 6.7e-16 (r01). |
| L4 | lift rows are relations | lift_n = primitive(b ш w), b ∈ B_n, w ∈ admissible words of weight 5−n. B1 = weight1.pickle: the branch rows 2(H1−H7)−3(H2−H6), H2−2H3+2H5−H6 (elementary: H(j)−H(8−j) = iπ(4−j)/4) and two Li_1 distribution rows; B2, B3 = A3c weight-2/3 selections (std + lifts); B4 = A3h5's T0[S] (std 4810, lift1 223, lift2old 7, lift3old 2) = Theorem W4R's inputs. b and w admissible ⇒ every word of b ш w is admissible ⇒ H(b ш w) = H(b)H(w) = 0 (shuffle product of absolutely convergent iterated integrals). "old" only means the bases are the earlier A3c/A3h5 selections. Own check: all 121630 lift rows regenerated with my own shuffle from B_n, identical label by label (r02a); all base and lift rows vanish numerically (≤ 5.4e-16 relative). |
| L5 | u = 7 image of a literal lift row is a relation | Rows have integer coefficients and live on admissible words; on [0,1] conj(z^j dt/(1−z^j t)) = z^{−j} dt/(1−z^{−j} t), so conj H(w) = H(7.w) and (7.r)·V = conj(r·V) = 0. Numerically: all 7-images of all lift families vanish (≤ 5.0e-16). |
| L6 | e ш W is a relation | e = (2) − (6), value iπ/2. W = Σ t_k G_k exactly (t05.log, 0 mismatches on 5184 words, 74/74), G_k = T0[S] rows and weight-4 lifts of 13 A3h4 frame rows. Each G_k is a proved relation (value 0 is used, not just span membership). Letters 2, 6 are admissible at both ends and G_k is on admissible words, so every word of e ш G_k is admissible and H(e ш G_k) = (iπ/2)·0 = 0. So e ш W = Σ t_k (e ш G_k) is a relation. **Weak point:** W itself (`wvecs.pkl`) was deleted, so the identity "the W of t04 is the W of t05" rests on the two logs of one run; I could not check it on stored data. My mod-p test (§3) replaces it by an independent existence check. |
| L7 | frame rows F_j (66 weight-5 rows, 15 distinct frame234 lifts with nonzero a_j, bases (2,3), (3,6)) and the 13 lower-weight rows | Function-level double shuffle of absolutely convergent multiple polylog series at real points of other Möbius frames, converted by `frame_ds.convert`. `candidates(n, bs, arity)` only yields factor tuples whose arguments come from bs = {1/t : t frame point, |t| > 1} (so real, |b| < 1), with |b| strictly decreasing along each factor and every vertex product (product of one cumulative b from each factor) again in bs; so every series in the stuffle product, and the product itself, is absolutely convergent. t06 asserts the stored factors are such a candidate and regenerates the row. "max |z| = 0.83" is the largest |ratio argument| x_i = b_i/b_{i−1} over the factors; it is a by-product, the convergence guarantee is the candidates() membership. Own check (r04): all 66 + 13 rows regenerated from provenance with my own domain test (real, 0<|x|<1, cumulative and vertex products = 1/t with |t|>1); 0 failures, max |x| = 0.8284. The 37 open rows' conjugate provenance has max |σ(x)| ≥ 1.707 (out of domain, as W5G1 said). The frame234 lifts equal primitive(base ш w) with my shuffle; lifts of a relation by an admissible word are relations (as L4). All 288 frame-w5 rows and all A3h4 rows vanish numerically (≤ 1.3e-16). Conversion code correctness = hypothesis (b). |
| L8 | t01–t03 (sectors, normal forms, rank-90 solve) | Used only to find a_j and the literal coefficients. Nothing they output is trusted by t04's final check except as numbers; an error there would make t04's check fail, not pass. Confirmed by reading t04: the sector data enter `target_coefs` only. |
| L9 | (H) ⇒ X·V = 0 ⇒ W_5 = N_5 + R_5 | X = Σ (W5C2 multipliers) u.f_i is W5G1 s05's definition; W5G1 s05 checked D = J + e ш Z + X exactly (0 mismatches, both components, in s05.log); W5G1R re-derived X and checked the split mod two primes, not exactly. J consists of literal rows and u = 7 conjugates; Z exactly in the proved weight-4 span (s10). |

## 2. What the result is conditional on

1. (a) The cited theorems: regularised (extended) double shuffle and distribution relations for level-8 multiple polylogarithm values; the shuffle product of convergent iterated integrals; Li_s(x) = (−1)^d G(…;1) and the stuffle product of absolutely convergent series. **And** the correctness of the code that turns them into rows (`relations.py`: `ds`, `distribution`, `reg`; and `cmzv.sh/mul`), exactly as in Theorem W4R. (My numerics check every row to ~1e-16, which is evidence, not proof.)
2. (b) Correctness of the frame conversion code (`frame_ds.identity/convert`, `cmzv.path/regularize`, Q(√2) arithmetic): now for 66 weight-5 and 13 weight-2..4 frame rows (W4R's caveat (iii)).
3. (c) W5C1's exact vector D (conversion of W_5, N_5, R_5 to words, including the path-reversal duality used for N_5).
4. (d) W5G1's exact identity D = J + e ш Z + X and W5C2's multipliers defining X. This identity was computed exactly by W5G1 (s05: 0 mismatches); the independent re-check (W5G1R) was modulo two primes only. So "exact" is justified as a statement about W5G1's computation; it is *independently* confirmed only modulo two primes.
5. The integrity of the W5G2 run itself: t04 and t05 exact checks are known from their logs; W was deleted, so the t04→t05 hand-over of W is not re-checkable without the 5-hour rerun. My independent check covers the mathematical content (existence of a certificate) mod two primes, not the exact certificate.
No hypothesis about Galois images (u = 3, 5) remains.

## 3. What I recomputed independently

- **Numerics (r01, own evaluator):** cumulative Gauss–Legendre on dyadic panels toward 0 and 1, all 52487 admissible words of weight ≤ 5 in 1 s; ζ(5) error 4e-16; H(1) and e = iπ/2 match mpmath. Weight-1 counterexample: 2(H1−H7)−3(H2−H6) = 0 (5.6e-16), its u = 3 image = 2πi (to 1e-15), u = 5 image = −2πi, u = 7 image 0. All 94762 std rows, all 121630 lift rows and their u = 7 images, all 288 frame-w5 rows, all A3h4 frame rows, B1–B4: max relative residual 6.8e-16. **All 74 targets u.f_i, built with my own action: max relative residual 1.5e-16.** Controls: letters 1↔2 swapped in u3.f4 gives 6.3e-3; u = 3, 5 images of the 10368 lift1 rows on branch base rows have |value| in [0.060, 8.33] (W5G1R saw [0.060, 6.8] on the 389 used ones).
- **Lifts (r02a):** all lift rows regenerated with my shuffle from B_n: lift1 20736, lift2old 31104, lift3old 34496, lift4old 35294; 0 mismatches against W5C1's rows; the 15 frame234 lifts with nonzero a_j: 0 mismatches.
- **Frame rows (r04):** see L7.
- **Weight 4 (r03):** G = T0[S] (5042) + weight-4 lifts of the 13 named A3h4 rows (94) = 5136 rows; rank 5070 mod 4194301 and mod 4194287 (t05 got 5070 with all 8686 generators, so the other A3h4 lifts add nothing); all 5136 complex conjugates 7.G lie in span(G) mod both primes; all G vanish numerically (5.3e-16). The frame rows t05 reports as used are exactly these 13.
- **W is not stored** (`cert_summary.pkl` holds only a_j and sizes; `wvecs.pkl` deleted). So I could not test "W ∈ proved weight-4 span" on W5G2's own W. Instead I tested the composite statement directly:
- **Span test (r02, all 74 targets, both sectors of u = 7, two primes):**
  x_t = u.f_i − Σ a_j F_j (a_j from `cert_summary.pkl`, F_j from data / my regenerated lifts), and the claim
  x_t ∈ span(L ∪ 7.L) + e ш span(G ∪ 7.G), L = std ∪ all lift rows (216392 rows).
  Method: L ∪ 7.L is stable under u = 7, so I project to the ± eigenspaces (dims 23382 / 23274); form n+64 random sparse combinations of the projected rows (their span is inside the true span, so a found membership is sound; randomness can only cause false negatives); own blocked dense elimination mod p carrying the test rows; then decide membership in L + e ш G on the residues (1200 random combinations of e ш G and, separately, of e ш words(4) = I).
  Results:

  | prime | sector | rank of span(L) | codim | rank of e ш G mod L | rank of I mod L | x_t ∈ L | x_t ∈ L + e ш G |
  |---|---|---|---|---|---|---|---|
  | 4194301 | + | 23113 | 269 | 7 | 72 | 0/74 | **74/74** |
  | 4194301 | − | 22993 | 281 | 21 | 70 | 0/74 | **74/74** |
  | 4194287 | + | 23113 | 269 | 7 | 72 | 0/74 | **74/74** |
  | 4194287 | − | 22993 | 281 | 21 | 70 | 0/74 | **74/74** |

  Second-prime rows checked by me against r02_4194287_{p,m}.log (finished 2026-11-23 10:00). Ranks and residue ranks are the same at both primes, so the random combinations reached the full span. The control pattern is also the same. In the − sector, u3.f4 itself (no frame terms subtracted) lies in L: its − projection does. That is harmless, because membership needs both sectors, and it fails in the + sector.

  Controls (each sector): a random vector and x(u3.f4) + one unit vector are outside L + e ш G and outside L + I; e ш (1,−1,−1,1) is in L + I but not in L + e ш G; u3.f4 and u5.f113 themselves (no frame terms subtracted) are outside L + e ш G in the + sector. So the test separates the cases, the frame terms are needed, and the e ш (proved weight-4) part is needed (no target lies in L alone).
  Reused project code in r02: none for algebra (own action, shuffle, projection, elimination); data reused: std rows (W5C1 `rels_std.pkl`, generated by `make_relations`), base rows B1–B4, frame rows, a_j. In r04 I reused `frame_ds.identity/convert/primitive` and `cmzv.POINTS/cross` (the conversion code under test is hypothesis (b)).

## 4. Minor remarks on result.md

- "23 weight-5 lifts … in the basis of 90": only 15 of them have a nonzero a_j in some certificate (bases (2,3) and (3,6)); the "66 frame-w5 + 13 A3h4 rows" count is right.
- Hypothesis (a) should name the row-generating code (`relations.py`) as W4R does; it is part of the standing caveat.
- Hypothesis (d): the split is exact in W5G1's computation; "the referee checked it modulo two primes" is accurate and should stay.
- The paper may say "proved, computer-assisted"; it should not imply the exact certificate was independently recomputed.

## 5. Paper text (paper/ieee/main.tex; each OLD string occurs once)

**(a) OLD**
```
For the other $37$ the images move the arguments outside the region of absolute convergence, and we have not proved that these $74$ images are relations; they vanish numerically to about $10^{-96}$. So $W_5-N_5-R_5$ is, exactly, a combination of these $74$ images (refereed once). The closed form of $W_5$ with three extra arguments is therefore computer-assisted and conditional on them (in our notes), the form with fourteen arguments is numerical, and $M_4$ is the first squircle coefficient that needs depth two (at the level of symbols).
```
**(a) NEW**
```
For the other $37$ the images move the arguments outside the region of absolute convergence, so the frame argument does not cover these $74$ images. We certified each of them directly. Each image is an exact rational combination of weight-five frame relations proved as above, of standard relations, of lifts, of complex conjugates of lifts, and of the shuffle product of $H(2)-H(6)$ with an exact rational combination of proved weight-four relations (those of the proof of Theorem~\ref{thm:W4R} and lifts of $13$ frame relations of weights two to four). That product is a relation, as above. We checked every identity in exact rational arithmetic on all $46656$ words, and the weight-four parts on all $5184$ words. No Galois argument is used. So the $74$ images are relations, and the closed form of $W_5$ with three extra arguments is proved, computer-assisted (in our notes). It rests on the same assumptions as Theorem~\ref{thm:W4R}: the cited double-shuffle and distribution theorems and the correctness of our code for the relations and for the frame conversion, which now produces $66$ weight-five and $13$ lower-weight frame relations. It also rests on the exact conversion of $W_5$, $N_5$ and $R_5$ to words. An independent re-check (refereed) regenerated the frame relations and lifts, found all $74$ images numerically zero, and confirmed their membership in the span of proved relations modulo two primes; it did not repeat the exact computation. The form with fourteen arguments is numerical, and $M_4$ is the first squircle coefficient that needs depth two (at the level of symbols).
```

**(b) OLD**
```
$M_4=W_5+W_4''+W_3^{(4)}+T$; $W_5$ needs depth two (symbol level) (Sec.~\ref{sec:square}-E) & $W_5$: Computer-assisted (three extras; conditional on $74$ Galois images of frame relations, checked numerically), Numerical (fourteen); $W_4''$: Computer-assisted; depth obstruction: Proved (symbol level) & 1596; exact certificate; exact symbol and cobracket \\
```
**(b) NEW**
```
$M_4=W_5+W_4''+W_3^{(4)}+T$; $W_5$ needs depth two (symbol level) (Sec.~\ref{sec:square}-E) & $W_5$: Computer-assisted (three extras), Numerical (fourteen); $W_4''$: Computer-assisted; depth obstruction: Proved (symbol level) & 1596; exact certificates; exact symbol and cobracket \\
```

**(c) OLD**
```
\item \emph{The next coefficient.} Prove the five weight-three and 29 weight-four level-8 relations (not needed for $M_3$) and prove that the images of $37$ weight-five frame relations under $\sqrt2\mapsto-\sqrt2$ are relations (the certificate for $W_5$ checks them only numerically), and extend the certificate to the form with fourteen arguments.
```
**(c) NEW**
```
\item \emph{The next coefficient.} Prove the five weight-three and 29 weight-four level-8 relations (not needed for $M_3$), and extend the certificate for $W_5$ to the form with fourteen arguments.
```

## Files

`rc.py`, `relim.py`, `r01_numerics.py` (r01.log), `r02a_lifts.py` (r02a.log), `r02_span.py` (r02_<p>_<p|m>.log/.pkl),
`r03_w4.py` (r03.log), `r04_frames.py` (r04.log). `H_values.pkl` and `my_lifts.pkl` were deleted after use; r01 and r02a regenerate them in seconds (r02 and r03 need them).
Rerun: r01 (10 s), r02a (15 s), r03 (2 min), r04 (15 s), r02 per prime and sector about 48 min, 5 GB.
