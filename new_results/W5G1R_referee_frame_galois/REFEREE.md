# W5G1R: referee report on W5G1 (Galois images in the W_5 certificate)

Scope: `W5G1_frame_galois/result.md` and its scripts, read only. All checks below are my own code in this folder
(`rcommon.py`, `r01`–`r08`). I use my own Galois action, my own shuffle, and my own double-precision evaluator
(graded panels with Chebyshev cumulative integration, ~1e-13). I did not re-run W5G1's scripts; where I cite their
logs I say so. Python: project `.venv`, python-flint 0.9.0.

Conventions (from `A3c_claude_elementary/Q/cmzv.py`): letter -1 = dt/t, letter j = z^j dt/(1 - z^j t), z = e^{2 pi i/8};
H(a1..an) = iterated integral over 1 > t1 > ... > tn > 0. Admissible: a1 != 0, an != -1. So H(j) = -log(1 - z^j).
The group acts by j -> u j mod 8, fixing 0 and -1. W5C2's g = 0,1,2,3 are u = 1,7,3,5 (`w5c2_fold.py`: tmap = -x, smap = 3x).

## Verdicts

| claim | verdict |
|---|---|
| 1 counterexample, W5C2 step 6 false | CONFIRMED after one small repair (value range) |
| 2 B, u = 7 | CONFIRMED |
| 3 C, std rows | CONFIRMED |
| 4 C, lifts, e sh Z | CONFIRMED |
| 5 A, frame rows | CONFIRMED with wording repair ("proved given the conversion code", same as W4R) |
| 6 remaining gap = X | CONFIRMED (one notation slip in W5C2, harmless) |
| 7 W4R and W4'' unaffected | CONFIRMED |
| 8 paper sentence | CONFIRMED AFTER REPAIRS (wording below) |

## 1. The value argument is false

By hand: for j = 1..7, H(j) - H(8-j) = log(-z^{-j}) = i pi (4-j)/4 (principal branch, argument in (-pi, pi)).
So 2(H1 - H7) - 3(H2 - H6) = 2(3 i pi/4) - 3(i pi/2) = 0. Its u = 3 image is 2(H3 - H5) + 3(H2 - H6) = i pi/2 + 3 i pi/2 = 2 pi i.
The u = 5 image has value -2 pi i. Numerically (`r01.log`, mpmath, 40 digits) the same.
This relation is row 0 of the weight-1 base used by `lift1` (`weight1.pickle`, selected rows). Row 1 is the j = 3 branch row.
Exactly: u.b0 = rho + 4e (u = 3) and rho - 4e (u = 5); u.b1 = rho - 4e (u = 3) and rho + 4e (u = 5), with rho in span(B1),
e = H(2) - H(6). The two distribution rows of B1 have e-part 0.

Usage, read directly from W5C2's four certificates (`r02.log`): 389 distinct `lift1` rows on base rows 0/1 are used.
Their u = 3 and u = 5 images enter the unfolded identity. With my own evaluator, word by word, over all 389 rows (`r06.log`):
|u.r.V| for u = 3, 5 lies in [0.060, 6.8]; for u = 1, 7 it is below 1e-14. My values reproduce W5G1's five printed ones
(0.0756, 0.156, 0.693, 0.590, 6.22). **Repair:** W5G1's range "0.08 to 6.2" comes from those 5 rows only; over all 389 it is 0.06 to 6.8.

So W5C2's step 6 ("(g.r).V = g(r.V) = 0") is false for u = 3, 5. The error is in the claim, not only in its proof.

## 2. B (u = 7)

Every row of every family is on admissible words, so each H(w) is a convergent integral on the real segment.
With t real, conj(z^j dt/(1 - z^j t)) = z^{-j} dt/(1 - z^{-j} t), and dt/t is real. So conj H(w) = H(7.w), and (7.r).V = conj(r.V) = 0.
Regularisation is used only to produce rows (`reg`, `regularize` in cmzv; `ds`, `distribution` in relations.py). Each final row
lives on admissible words and is an identity among convergent values. Nothing needs to commute with conjugation.
Also, all 94762 std rows are literally closed under u = 7 (`r03.log`), and the u = 7 image of each of the 38 used frame rows
is literally a frame-w5 row (`r05.log`).

## 3. C for std rows

The weight-5 std family is `make_relations(5, {})`: regularised double shuffle (92087) and distribution (2675).
It has no branch rows, no separate "regularisation rows" and no lifts. Why it is closed:
- `to_series` uses differences of consecutive letters; these scale by u. The stuffle merge adds the second entries mod 8,
  which commutes with scaling. `from_series` takes cumulative sums. Shuffle commutes with any letter map.
- `reg` singles out only the letters 0 and -1, which u fixes. So reg(u.w) = u.reg(w), and the base-point convention
  H(0) = H(-1) = 0 is invariant.
- The pairs (u, v) range over all words (first factor non-admissible allowed) and ds is symmetric, so the set of generated
  rows is closed.
- Distribution d in {2, 4, 8}: words with letters in d Z/8 go to such words, the lifted letter sets
  {j/d + k 8/d} go to themselves since u is a unit mod d, and the set C of nonzero multiples of 8/d is u-stable.
Computation (`r03.log`): all 94762 rows, images under u = 3, 5, 7 are literally std rows (up to the primitive scaling). 0 failures.
Note: the lower-weight bases are not closed (the branch rows are not); that is handled in claim 4, not here.

## 4. C for lifts, e sh Z

Logic. A lift row is r = kappa (b sh w) with b a proved base row. Then u.r = kappa (u.b sh u.w) and u.b = rho + e sh y exactly,
rho in the proved base span. So u.r = kappa (rho sh u.w) + e sh (kappa y sh u.w). The first part is a shuffle of a proved relation
with a word, hence a relation. The second parts add up to e sh Z. If Z = sum t_i G_i with G_i proved weight-4 relations,
then e sh Z = sum t_i (e sh G_i). Every word in e sh G_i is admissible (e has letters 2, 6; G_i is on admissible words),
so H(e sh G_i) = H(e) H(G_i) by the shuffle product of convergent iterated integrals on [0,1] (no regularisation),
= (i pi/2)(0) = 0. So e sh Z is a relation. The value map is compatible with "shuffle with e" here precisely because
all words are admissible.

Checks:
- The 56 stored decompositions u.b = sum t_i b_i + e sh y (`lift_decomp.pkl`) hold exactly with my shuffle; the only nonzero
  e-parts are the four weight-1 ones above (`r03.log`).
- The bases are proved: B1, B2, B3 are regenerated by `make_relations`; B4 = T0[S] (std 4810, lift1 223, lift2old 7,
  lift3old 2) matches regeneration from its labels, 0 mismatches (`r03.log`).
- Z recomputed from scratch from the certificates' exact y, the lift1 rows and the hand decomposition: supports 662 / 666,
  as W5G1 (`r04.log`).
- With W5G1's multipliers (`Z_certificate.pkl`) and generators rebuilt by me (B4 from r03; the A3h4 frame rows regenerated
  from provenance through frame_ds with my own domain checks, max |z| = 0.828): sum t_i G_i == Z exactly on all 5184 words,
  both components, 0 mismatches; also equal on random projections mod 2^61 - 1 and mod 1e9+7.
- The 19 frame lifts in Z are exactly the 19 weight-2/3 frame lifts of Theorem W4R's certificate (frame rows (2,3), (3,3),
  (3,5), (3,6) with the same words). So Z uses nothing beyond W4R's inputs.

## 5. A: the 38 frame rows

I regenerated all 38 used rows, plus rows 160 and 286 (the literal u = 3, 5 images of row 80), from stored provenance with
`frame_ds.identity/convert` and my own domain checks (arguments real, 0 < |z| < 1, cumulative and vertex products equal
1/t for frame points with |t| > 1). All 40 pass and equal the stored rows; max |z| = 0.8284 = 2 sqrt2 - 2 (`r07.log`).

What "proved" rests on, exactly:
(a) theorems: product of absolutely convergent multiple series = stuffle sum; shuffle product of convergent iterated integrals;
Li_s(z) = (-1)^d G(...; 1);
(b) correctness of the conversion code: Möbius change of variables, Cayley pullback dlog((x-p)/(x-q)), homotopy of the image
arc to the radial path through the disk centre, piecewise shuffle regularisation in `cmzv.path` whose total is a convergent word,
and exact Q(sqrt2) arithmetic. This is audited (`A3h4_frame_ds/PLANNER_AUDIT.md`) and validated numerically (W5C1: 2.8e-113),
but not formally verified.
So the fair status is "proved, given the correctness of the conversion code", the same as the 28 frame lifts of Theorem W4R.
Bare "proved" is too strong unless the paper's standing caveat (iii) of Theorem W4R is stated to cover it.
W5G1 is right that W5C2's verifier did not regenerate frame rows; s06 and r07 now do.

## 6. The remaining gap

- W5C2's unfolded identity: my mod-p check (two primes, random projections, own code) passes with c = y/(4 lam) for both
  components. It fails with c = y lam/4. W5C2's result.md writes "y_i lam_i/4"; its code (and s05) uses y/(4 lam).
  Lambda takes the values 1, -1, 2, 4, so the text is a notation slip. The certificate itself is fine.
- X recomputed by me: the 74 non-literal pairs are rows {4, 5, 7, ..., 113} (the 37 listed) x {3, 5}; row 80's images are
  literally rows 160 and 286. Supports 26286 / 26264, as W5G1.
- D - X - e sh Z equals, on random projections mod both primes, the sum of proved vectors only: literal rows, u = 7 images,
  std images (literal std rows), lift images rebuilt as kappa sum t_i (b_i sh u.w) from the proved bases, and the two literal
  frame images (`r05.log`). Controls: X.pi and (e sh Z).pi are nonzero mod p (`r05b.log`), so the test sees both.
- I did not repeat W5G1's exact full-vector reassembly (s05, 22 min). With s05's exact check, W5C2's own exact step_05,
  and my independent mod-p checks, I accept D = J + e sh Z + X.
- Hidden ingredients in J: none new. J uses only the premises already in Theorem W4R's list (std theorems, shuffle,
  conversion code) plus the target vector D itself (W5C1's conversion of W5, N5, R5; N5 uses the path-reversal duality
  for Li_{3,2}(1/sqrt2, 1)). Those are outside the Galois question.
- Numerics: with my double-precision evaluator, all 74 u.f vanish to 3.4e-12 relative and X to 2e-13 relative; a control
  (letters permuted) gives 0.16 (`r08.log`). W5C2's planner checked all 38 frame rows x 4 g at about 100 digits: max 1.6e-98
  (`W5C2_exact_certificate/planner/planner_verify_w5c2.out`). That supports the paper's "about 10^{-96}". W5G1's s07 (1e-61)
  is lower precision, not a contradiction.
- The mod-p C2 membership (s09) I did not re-run. It is one prime, in sector quotients, and membership in proved + e sh (anything)
  does not imply being a relation. It is consistent with the motivic heuristic, not evidence of much more.

So: W5 = N5 + R5 + X.V holds as a theorem (same standing caveats as W4R), and W5 = N5 + R5 is conditional on X.V = 0,
i.e. on the 74 Galois images (37 frame rows x {3, 5}) being relations, or at least on their weighted sum X vanishing.

## 7. Theorem W4R and W4''

- `A3h5_exact_certificate/verify_D.py` rebuilds every row: std from `make_relations(4, {})`; lift1/lift2old/lift3old as
  primitive(b sh w) with b from weight{1,2,3}.pickle (asserted to be regenerated rows); frame rows from provenance through
  candidates/identity/convert with domain asserts. There is no letter map, no sector, no fold, no conjugation. The lift1 base
  rows are the literal branch rows (j = 1 and j = 3, each generated directly by its own branch formula) and two distribution rows.
  Families used: std 4810, lift1 223, lift2old 7, frame 28 (10 + 9 + 9 lifts of 13 distinct frame rows).
- W4'' is `M4b_w4pp_exact` (certificate `W4pp_certificate.json`, verifier `verify_W4pp.py`) and the nicer form
  `M4b2_w4pp_nice` (`verify_W4pp_nice.py`). Both verifiers are the same template (diff: certificate path and log text only).
  Families: atom 29, std 4810, lift1 223, lift2old 7, frame lifts 28 (the same 10 + 9 + 9). No sectors, no Galois images.
  grep for conj/sigma/letter maps in the verifiers, frame_ds.py, relations.py, cmzv.py: none used.
- So the u = 3, 5 error does not apply to W4R or W4''. Both still carry the same caveat (iii): correctness of the conversion code.

## 8. Paper wording

Problems with W5G1's proposed sentence:
- "their total is i pi times proved weight-four relations, in shuffle product" is muddled. The defect is the shuffle product
  of H(z^2) - H(z^6) (value i pi/2) with an exact combination of proved weight-four relations; that product is a relation.
- "found them in the span of proved relations and of i pi times weight-four words modulo a prime" should say "one prime" and
  should not suggest more than it shows (or be dropped).
- The earlier sentence "In each sector the vector is an exact rational combination of standard relations ... and ... 36
  weight-five function-level double-shuffle relations" hides the Galois images. It should say so. 38 distinct frame relations
  are used, 36 in each of two sectors.
- "proved in the same way as at weight four" is fine only if caveat (iii) of Theorem W4R is understood to apply.

Suggested replacement, from "The Galois group of ..." to the end of that paragraph's W_5 part:

"The Galois group of $\mathbb Q(e^{2\pi i/8})$ over $\mathbb Q$ splits the problem into four sectors. Unfolded, the
certificate writes $W_5-N_5-R_5$ as an exact rational combination of relations and of their images under the Galois group,
which acts on the letters by $j\mapsto uj$, $u\in\{3,5,7\}$. The relations are standard ones (regularized double shuffle and
distribution), lifts of lower-weight relations, and $38$ weight-five function-level double-shuffle relations in other M\"obius
frames ($36$ in each of two sectors). The multipliers come from an exact $p$-adic linear solve (heights up to about $2080$
digits), and the check is exact rational arithmetic on all $46656$ words. Images under $u=7$ (complex conjugation) are
relations. Images under $u=3,5$ (both send $\sqrt2$ to $-\sqrt2$) need not be: the relation
$2(H(1)-H(7))=3(H(2)-H(6))$ is sent to a combination with value $2\pi i$. The images of the standard relations are again
standard relations. The images of the lifts differ from proved relations by the shuffle product of $H(2)-H(6)$, whose value
is $i\pi/2$, with an exact rational combination of proved weight-four relations (those of the proof of
Theorem~\ref{thm:W4R}); this product is again a relation. The $38$ frame relations themselves are proved as in
Theorem~\ref{thm:W4R}, with the same conversion code and the same convergence conditions. For one of them its images are
again frame relations. For the other $37$, the images move the arguments outside the region of absolute convergence, and we
have not proved that these $74$ images are relations; they vanish numerically to about $10^{-96}$. So
$W_5-N_5-R_5$ equals, exactly, a combination of these $74$ images, and the closed form of $W_5$ with three extra arguments
is computer-assisted and conditional on them (in our notes)."

Status labels a paper may print:
- Proved (computer-assisted, with Theorem W4R's caveat (iii)): Theorem W4R; W4''; the 38 weight-five frame relations;
  "W_5 - N_5 - R_5 is an explicit rational combination of 74 Galois images of frame relations".
- Conditional: W_5 = N_5 + R_5 (three extra arguments), on those 74 images being relations (or on their combination X vanishing).
- Numerical: each of the 74 images vanishes to about 1e-96; the form with fourteen arguments.
- Heuristic only: that Galois images of relations hold modulo shuffle with 2 pi i (motivic argument, not checked or cited).
Table: "Galois images of 37 frame relations checked numerically" is right. Open problem: "prove that the images of 37
weight-five frame relations under the automorphisms sending $\sqrt2$ to $-\sqrt2$ are relations" is right.

## Files

`rcommon.py` (action, shuffle, evaluator), `r01_weight1.py` (claim 1 by hand and numerically), `r02_extract.py` (reads W5C2
certificates), `r03_closure.py` (std closure, base provenance, decompositions), `r04_Z.py` (Z and its certificate),
`r05_identity_modp.py` + `r05b_controls.py` (unfolded identity, X, reassembly), `r06_lift_values.py` (389 x 4 values),
`r07_frames.py` (frame regeneration, conjugate arguments), `r08_numerics.py`. Logs r0N.log. Run in this folder, in order;
total under 5 minutes.
