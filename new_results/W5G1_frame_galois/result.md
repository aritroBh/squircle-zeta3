# W5G1: Galois conjugates in the W_5 certificate

## First: an error in the existing certificate

Step 6 of `W5C2_exact_certificate/result.md` (the "value argument" (g.r).V = r.(gV) = g(r.V) = 0) is false, not just unproved.
For u = 3, 5 the conjugate of a true relation need not be a relation. The smallest example is weight 1:
2(H1 - H7) - 3(H2 - H6) = 0 is a true relation (it is a row of the weight-1 base used by `lift1`), but its u = 3 image
2(H3 - H5) + 3(H2 - H6) has value 2*pi*i (s02, s03).
W5C2 uses the u = 3, 5 images of `lift1` rows built on that base row and on one other branch row. They are not relations.
Numerically their values are 0.08 to 6.2 (`s02.log`); 389 used rows x 2 values of u are of this kind (`s05.log`).
The planner's numerical check sampled 40 non-frame rows and happened to miss them.

So W5C2's unfolded identity D = sum (y lam/4) sum_g chi(g) g.r is a correct exact vector identity (s05 re-derives it,
0 mismatches), but some of its terms are not relations. As written, it does not prove W5 = N5 + R5.
The repair below proves everything except the conjugates of 37 frame rows.

## Status lines

- **A (frame rows themselves): PROVED GIVEN the audited conversion code**, the same status as weight four.
  `w5_frame_gen2.py` calls `frame_ds.frames/candidates/identity/convert` unchanged, the same code that produced the 28 lifts in Theorem W4R.
  s06 regenerates all 38 used frame-w5 rows from their stored provenance and re-checks the domain conditions: real arguments,
  |z| < 1 (max |z| = 0.828), and every cumulative and vertex product a frame point with |t| > 1. All 38 equal the stored rows.
  One difference from weight four: W5C2's verifier reads `validated_w5.pkl` and does not regenerate the rows from provenance; s06 now does that.
- **B (u = 7): PROVED.** Every vector lives on admissible words (first letter not 0, last letter not -1), so each H(w) is a
  convergent iterated integral over the real segment [0,1]. The forms are dt/t and zeta^j dt/(1 - zeta^j t). With t real,
  conj(zeta^j dt/(1 - zeta^j t)) = zeta^{-j} dt/(1 - zeta^{-j} t), so conj H(w) = H(7.w).
  For a rational row r this gives (7.r).V = conj(r.V) = 0. Regularisation only enters the generation of rows. The rows
  themselves are identities among convergent values, so no convention has to commute with conjugation.
  This holds for every row family. In addition, the frame-w5 set is literally closed under u = 7, because x -> -x preserves the real point set and |z|.
- **C (u = 3, 5 on the frame rows): OPEN** (rigorously), with exact partial results and mod-p evidence.
  - Lifts and std (not frames): **PROVED, exactly.** The u = 3, 5 images of std rows are literally std rows.
    For a lift row b sh w, s04 shows exactly that u.b = rho + e sh y, where rho is in the proved base span, e = (2) - (6) and e.V = i*pi/2.
    A nonzero e-part occurs only for the two weight-1 branch rows.
    s05 collects all defects into e sh Z, with Z a weight-4 vector. s07/s10 show exactly that Z lies in the span of
    5056 proved weight-4 rows: A3h5's T0 rows plus 19 lifts of the A3h4 frame rows (2,3), (3,3), (3,5), (3,6).
    Z is not in T0's span alone. Hence e sh Z is a sum of shuffle products of proved relations, and is proved.
  - C1 (structural): u = 5 is sigma (sqrt2 -> -sqrt2) on the real points, because the Cayley exponents satisfy x -> 5x. u = 3 is sigma followed by x -> -x.
    The conjugate provenance is the same function-level identity in the conjugate frame at the conjugate points. The conjugate frame is always one of the 336 frames.
    - Row 80 (u = 3 and u = 5): **literal.** The conjugate is in the domain, and regenerating it gives exactly u.row.
    - The other 37 rows x {3, 5} fail the domain check, because conjugate arguments such as -(2+2*sqrt2) or 2+sqrt2 have |z| > 1 (`s06.log`).
  - C2 (span, mod p = 16777213, one prime; NOT exact).
    - Let Et = span(std, all lift rows, I), with I = e sh Q^{words(4)}. Et is Galois-stable, and its sector quotients are computed in s08: nullities 98/99/107/104, total 408.
    - In that quotient, every one of the 74 vectors u.f, and the aggregate X, lies in Et + span(frame-w5 rows, frame-234 lifts, u = 7 images) (s09).
    - The frame-w5 rows alone do not suffice: the rank goes 67 -> 72. The frame-234 lifts are needed.
    - Caveat: I is not made of relations. Membership in Et means u.f = (proved) + e sh Y + (frame terms), and Y must still be shown to be a proved weight-4 relation, as Z was.
    - Not done: an exact rational certificate (Dixon solves in the four sectors, plus a weight-4 check on Y).
  - Numerics (s07, 60 digits): all 76 u.f.V are about 1e-61 (precision-limited); X.V/max|X| is about 1e-59; Z.V is about 1e-58.
- **D (Theorem W4R / A3h5): PROVED, no gap.** A3h5's final identity is a direct rational combination of 5068 literally generated rows on all 5184 words:
  std 4810, lift1 223, lift2old 7, frame 28. There are no sectors and no conjugates.
  `verify_D.py` regenerates every row, including frame rows from provenance with domain asserts, and contains no Galois map.
  The lift1 rows there are lifts of true weight-1 relations, not their conjugates. Not checked here: W4'' ("a certificate of the kind of Theorem W4R"); if it used sectors, it should get the same audit.
- **E (rebuild): PARTIAL.**
  - Exactly, on all 46656 words and for both components: D = J + e sh Z + X (s05, 0 mismatches), with Z = sum t_i G_i exactly (s10).
  - J is a combination of 77492 literal rows and u = 7 conjugates of literal rows; the G_i are proved weight-4 rows.
  - So D - X is an exact combination of justified relations. X is the u = 3, 5 images of 37 frame rows with W5C2's multipliers (support 26286 / 26264).
  - Residual: X, not 0. It is unresolved exactly; it is in the proved span mod p, up to the I-caveat above.

## (row, u) pairs by route (38 used frame-w5 rows, indices into validated_w5.pkl)

- u = 1 and u = 7: all 38 rows, literal or by B.
- u = 3, 5, literal / C1: row 80.
- u = 3, 5, C2 (mod p only, not exact): rows 4, 5, 7, 9, 10, 12, 14, 15, 16, 17, 52, 53, 54, 56, 58, 59, 61, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 87, 91, 95, 103, 104, 107, 113.
- Exactly unresolved: those same 37 x 2.

## C3: what would suffice, and what the literature gives

- The naive statement "g in Gal(Q(zeta8)/Q) maps relations to relations" is false (weight-1 example above).
- What transport of structure gives: the 8th roots of unity form a Q-rational set, so sigma acts on the motivic fundamental groupoid of P^1 - {0, mu8, inf} with tangential base points 0 and 1.
  - It sends the de Rham letters omega_a to omega_{ua}.
  - It acts on de Rham periods, the canonical-fibre-functor periods, i.e. motivic periods modulo 2*pi*i.
  - So if r holds for motivic iterated integrals (Deligne–Goncharov), then u.r holds modulo the ideal generated by 2*pi*i.
- That is exactly the e sh (.) defect seen here (unverified heuristic, my own; no citation located).
- What would suffice, labelled as a hypothesis: (H) for each of the 37 frame rows f, u.f - (proved) lies in e sh (proved weight-4 relations).
  - The s09 mod-p result is consistent with (H).
- On the literature, as far as I know (not checked against the texts in this job):
  - Racinet proves double shuffle and distribution for the numbers.
  - Deligne–Goncharov construct the motivic groupoid.
  - Glanois works with motivic iterated integrals at N = 8.
  - Souderès proves motivic double shuffle at N = 1.
  - I know of no published motivic version of function-level double shuffle at real algebraic points in Möbius frames, so (H) cannot be cited.

## Proposed replacement for the paper sentence

Replace "The proof uses the same ingredients as Theorem thm:W4R ... but did not prove either." with:

"The weight-five frame relations are proved in the same way as at weight four. The sector decomposition also uses the images of the relations under $\sqrt2\mapsto-\sqrt2$. For the standard relations and the lifts these are handled exactly: the images of lifts of the two weight-one branch relations are not relations themselves, but their total is $i\pi$ times proved weight-four relations, in shuffle product. What we have not proved is that the images of 37 of the 38 frame relations used are relations; their arguments leave the region of absolute convergence. We checked them numerically to about $10^{-96}$ and found them in the span of proved relations and of $i\pi$ times weight-four words modulo a prime, so the closed form of $W_5$ with three extra arguments is computer-assisted and conditional on them (in our notes)."

Also adjust:
- table line ~2044, "frame relations checked numerically", to "Galois images of 37 frame relations checked numerically";
- open problem line ~2200, "prove the weight-five frame relations, and their Galois conjugates", to "prove that the images of the 37 weight-five frame relations under $\sqrt2\mapsto-\sqrt2$ are relations". The frame relations themselves are proved.

## Files and rerun (python = ../../../.venv/bin/python, run in this folder, in order)

- s01_used_rows.py: used rows and y from W5C2's certificates.
- s02: numerical counterexample.
- s03: low-weight defects, mod p.
- s04: exact decompositions → lift_decomp.pkl.
- s05: exact split D = J + e sh Z + X, about 22 min → split.pkl (750 MB).
- s06: A plus C1.
- s07: Z membership plus numerics.
- s08: sector quotients, about 2.2 h → quot_*.npz.
- s09: C2 tests.
- s10: exact Z certificate → Z_certificate.pkl.

Logs are in s0N.log/.out. used_rows.pkl (357 MB) and split.pkl can be regenerated by s01/s05.

## Planner notes after the referee (2026-11-22)

Referee: `../W5G1R_referee_frame_galois/REFEREE.md` (independent, own Galois action, shuffle code and evaluator). All eight claims confirmed. Repairs:
1. The value range of the 389 false lift images is 0.060 to 6.8 over all rows (the "0.08 to 6.2" above is from 5 sampled rows).
2. W5C2's unfolded identity has coefficient y/(4λ), not yλ/4 (notation slip in W5C2's result.md; the certificate is fine).
3. "Proved" for the 38 frame rows means "proved given the conversion code", the same caveat as in Theorem W4R.
4. Count: 38 distinct frame relations, 36 in each of two sectors. Open: 74 images (37 rows × u = 3, 5); numerically zero to 1.6e−98.
5. The mod-p span membership is one prime only and proves nothing by itself.
Theorem W4R and both W₄″ certificates (M4b_w4pp_exact, M4b2_w4pp_nice) use no sectors and no Galois map: not affected.
The two large intermediates split.pkl (750 MB) and used_rows.pkl (357 MB) were deleted after the referee finished; s01 and s05 regenerate them (about 25 minutes).
