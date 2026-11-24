# Planner audit of the frame double-shuffle rows (A3h4): domains, paths, branches

Audited by reading `frame_ds.py` (Codex). It answers the top session's condition 1: for each frame row, the region where
the identity holds, the branch conventions, and why the specialised point is inside the region.

## The identity (`identity(factors)`)
For factors u = Li_{s}(z) given as series words, the identity is
  (product of the series, expanded by the stuffle) − (product of the iterated integrals, expanded by the shuffle) = 0,
written entirely as G-words via Li_{s1..sd}(z1..zd) = (−1)^d G(0^{s1−1}, 1/b1, …, 0^{sd−1}, 1/bd; 1), where
b_j = z1···zj (`gword`). `numerical_selfcheck` checks this convention and the identity on 25 instances, against direct
series, to 1e-66.
- Stuffle side: the product of absolutely convergent multiple series is the stuffle sum. Absolute convergence allows
  the rearrangement.
- Shuffle side: the product of convergent iterated integrals over the same path is the shuffle sum (Chen).
So each row is a theorem wherever all series converge absolutely and all iterated integrals converge on the path.

## The domain enforced by the code
- Frame: a real Möbius map, t = cross(x; p0, p1, pinf), with p0, p1, pinf among the 8 real configuration points
  {0, ±1, ±A, ±B, ∞} (`frames`). Its frame points ts are the t-images of the 8 points, all real.
- Admissible cumulative products: bs = {1/t : t a frame point, |t| > 1}, so every b ∈ bs is real with 0 < |b| < 1.
- `sequences`: for each factor, the cumulative products satisfy |b_d| < … < |b_1| < 1. So every argument
  z_j = b_j/b_{j−1} has |z_j| < 1 (and |z_1| = |b_1| < 1).
- `candidates`: every "lattice vertex", i.e. every product of cumulative b's across factors, is required to lie in bs.
  The stuffle merges arguments by multiplication, so every cumulative product in every stuffle term is real with
  modulus < 1, and every G-letter 1/b is a frame point with |1/b| > 1.
- Consequence: every Li-series involved converges absolutely; strict |·| < 1 means no (s1, z1) = (1, 1) divergence.
  Every G-word has its innermost letter 1/b_d ≠ 0 (convergent at t = 0) and its outermost letter ≠ 1 (convergent at
  t = 1). All non-zero letters are real and outside [−1, 1], so the straight path [0, 1] meets no singularity.
  No analytic continuation is used and no branch of log or polylog is ever chosen. The points lie in the open region of
  absolute convergence, not on its boundary.

## The conversion to level-8 words (`convert` → cmzv.forms_path / path)
- dt/(t − a) = dlog((x − g(a))/(x − pinf)), with g = t⁻¹. This is exact for a Möbius map with pole at pinf
  (`dlog(ts[a], pinf)`); the letter a = 0 gives g(0) = p0.
- The t-path [0, 1] maps to the arc of RP¹ from p0 to p1 that avoids pinf and every letter g(a) (none of the a lie in
  [0, 1]). Under the Cayley map u = (x − i)/(x + i) this arc lies on the unit circle and avoids all letters.
- All letters of these forms lie on the unit circle (dlog((x−p)/(x−q)) = dlog((u−u_p)/(u−u_q)), with u_∞ = 1). So the
  open disk contains no singularity. Pushing the arc into the disk gives a homotopy, avoiding every letter, to the
  engine's radial path u(p0) → 0 → u(p1). Hence the iterated integrals agree.
- Inside `path`, the radial pieces are shuffle-regularised individually. Their total is a convergent word, so the
  regularisation cancels. This engine was validated on the W4/R conversion (1e-50 to 1e-51).
- Where an iπ-type error could hide (point raised by the top session). The homotopy argument alone does not justify
  the regularised path composition inside `path`.
  - For frame rows, u = 0 is NOT a letter: every form is dlog((u − u_p)/(u − u_q)) with u_p, u_q on the unit circle. So
    the radial path through 0 is a genuine path there, with ordinary Chen composition.
  - The delicate places are the endpoints u(p0) and u(p1). After splitting at 0, the individual pieces can have
    letters equal to an endpoint, i.e. they diverge there, and they are shuffle-regularised with the engine's
    tangential-base-point convention. The divergences cancel only because the total word converges and the same
    convention is used in every piece.
  - That step therefore rests on (a) the engine's validated path convention (the A3c conversions of W4, R and their
    atoms agree to 1e-50 to 1e-51), and (b) the per-row numeric guard below. A missing iπ·(lower-weight) term would be
    O(1), not 1e-101.
- Sanity that no stray branch term entered: every frame row vanishes numerically, to 1e-101 (Codex: 44 / 92 / 184
  rows at weights 2 / 3 / 4), and the weight-3 rows also vanish to 1e-139 against the INDEPENDENT 150-digit values
  of A3f6 (planner_check_w3_numeric.py). A branch error would show up as an O(1) term such as 4πi·log(1+√2), the
  failure mode found in A3h2's Galois test.

## Lifts
Frame rows are lifted to weight 4 by shuffling with convergent words (`words(4−n)`, admissible). The product of
convergent iterated integrals from 0 to 1 is their shuffle (theorem).

## Verdict
Every frame row and frame lift is a specialisation of the function-level double shuffle (absolutely convergent
multiple series ↔ Chen iterated integrals), inside the open region of absolute convergence, followed by an exact
homotopy-invariant change of variables. No branch choices are involved. PROVED-status rows, subject only to the
correctness of the code that writes them down. Condition 2 below guards that.

## Verifier (the top session's condition 2)
`verify_D.py` rebuilds each frame row it uses from the function-level identity: it recomputes `identity(fs)` with
`convert(g, ts, ·)` from the recorded (frame, factors), asserts that fs ∈ candidates(...) (the domain checks above),
and compares the result with the stored row. T0 rows are rebuilt from relations.make_relations, and D from the
original W4/R sources. Remaining caveat: the verifier's correctness rests on frame_ds.py and cmzv.py, which this
audit read but did not re-implement.
