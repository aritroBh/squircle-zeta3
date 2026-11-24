**Overall verdict: PROOF COMPLETE, with the wording corrections below, for the tower law through depth 10 and the weight-seven law for every integer k >= 7.** Task E's finite non-adjacent reduction is correct; a broader all-depth tower equivalence needs additional premises, and its denominator-free prefactor prose is erroneous. These outlook issues do not affect the finite-depth theorem.

CX11 independent referee report. Admissible depths are 0 <= delta <= min(10,2d-4); off-diagonal levels are k >= d; stable diagonal levels are k >= 2p+1. The actual initial diagonal value at k=2p is retained separately. I read TOW5 `result.md` and `proof_depth10.md` fully, all requested Python sources, the inherited S-derivative derivation, CX10 and TOW4R reports, and WTS Lemmas 2–3. All writes were confined to this folder. Only executed checks are described as fresh passes.

1. **OK — R_5, its symbolic method, and the depth-10 seed.**

   `prove_S5.py` reran with exit 0: Stirling recurrence and anchors, symbolic rational identity, reduced denominator 3 Pi_5, degree-10 numerator, and all 24 supplied direct derivatives pass (`r5.log`). Its six finite-sum functions `ff`, `tf`, `tf1`, `normalized_summand`, `Pi`, and `total` are AST-identical to `prove_S34.py`, which is byte-identical to TOW2's source. Extending the finite bounds to l=5 requires only the new tail P_5. The 126-term sum keeps d and m symbolic; N_5 is computed from that sum, with no fit. Intermediate falling-factorial divisions must be cancelled before specializing to boundary keys. The final Pi factors are nonzero throughout the depth-10 domain d>=7.

   Independent `check_r5.py` derives P_5 by counting nonsingleton block multiplicities with sum (j-1)c_j=5, using falling(x,sum j c_j)/prod((j!)^c_j c_j!). It reproduces the stated polynomial exactly and checks Stirling values at n=6..35. Direct exact fifth derivatives of independently constructed S polynomials match at fresh (d,m)=(11,11),(11,7),(13,8),(14,7),(17,12), including an axis and a diagonal (`independent_r5.log`). The freshly regenerated N_5 equals the submitted N_5 symbolically (`shapes.log`).

   **Seed-blindness needs precision.** With lower depths fixed, adding a constant to a_10 leaves closing unchanged. Adding d^5/Pi_5 also preserves closing while changing the top by 120/Pi_5. Thus closing cannot independently certify the seed or identify the prescribed zigzag top. It is not blind to every rational seed error: adding d^6/Pi_5 changes closing by 720/Pi_5 (`seed.log`). The needed conclusion is valid: rho(10) rests on the proved R_5 -> sl_5 -> L(10) -> Bell derivation, not on closing alone.

2. **OK — a_9, a_10, tau, Bell, and fresh towers.**

   The original `verify_L.py` passed r=1..8 before interruption for expression growth. An equivalent rerun first cancelled the cumulants in Q(i)(d,m), then executed all ten original sqrt(2)-aware numerator comparisons: every r=1..10 passes, exit 0 (`L_reduced.log`). This checks equality to the proved L(r) construction, rather than agreement between two candidate caches.

   `check_initial.py` passes (`initial.log`): the defining tau expressions give tau1(9)=(-99-297i)/256, tau2(9)=(99-297i)/256, tau1(10)=1-99i/256, tau2(10)=1+99i/256, exactly; conjugation symmetry holds. The independent recurrence n a_n=sum_{r=1}^n r L_r a_{n-r} agrees with partition/Bell assembly at every depth 1..10. A separate cumulant recurrence reproduces sl_5. Both saved a_9/a_10 closed forms equal Bell symbolically (`details.log`). Since t'(i)=0, depth 9 needs only R_1..R_4, while depth 10 first needs R_5.

   Byte-copied TOWN1 `trunc_engine.py` and `driver.py` generated new exact Laurent series here. For d=7,8,11,13 and p in {d,d-1,ceil(d/2)}, all depths 0..10 satisfy rho(delta)=rho(0)a_delta: **132/132 equalities**. Fresh tower propagation at k=14 and 16 gives **20 further depth-9/10 equalities**, at keys (7,0),(4,3),(5,3),(4,4),(7,4), using stable ehat for (4,4). These checks regenerate series and towers; they do not read TOW5's saved tower cells. The off-diagonal factor is already in rho(0).

3. **OK — forced recursion, closing, and uniqueness.**

   A genuine rebuild began with only the depth<=8 cache; the submitted depth-9/10 cache was removed locally before execution. `forced910.py` proves F(9,5)=0 in the raw run (`forced.log`) and F(10,6)=0 in the equivalent exact-field rerun (`forced_reduced.log`, exit 0). The latter retains the genuinely rebuilt depth-9 family. Every rebuilt entry equals the submitted cache symbolically, all depth<=6 entries equal the inherited TOW3 family, and every top through depth 10 equals its zigzag formula (`details.log`).

   The recursion is ascending in Newton index s. Every equation below the top holds by construction; closing supplies the top equation, and higher indices vanish by support. Adjacency and E3 remain separate obligations. At fixed finite k, start at degree k with actual initial data and descend in degree using the triangular tower recursion. At a diagonal initial key retain actual rho; E3 permits its feed to the adjacent child to be replaced by stable extrapolation. This finite recursion has a unique solution, proving the constructed candidate equals the actual tower. No interpolation or induction from infinite degree is used.

4. **OK — all eleven adjacent identities, slow validation, and mutation.**

   All cases **(9,0..4),(10,0..5)** pass with exact zero numerator and exit 0 (`adjacent.log`, individual `adj_*.log`). The diagonal basis conversion is Ediag(delta,s)=sum_t(-1)^t ehat(delta,s+t). The integer-power reductions in the fold use i^(4p)=1 for integer p before rational continuation; they do not assert that equality for arbitrary complex p.

   `check_method.py` adapts the original TOW3 slow expressions only by extending the family/tau providers. At **(9,3)** and **(10,4)** the lhs and every individual feed term agree with the fast expressions (30 and 33 terms), and independent summation of the slow expressions in Q(i)(p) gives zero residual. Perturbing the actual family input **F(9,4) by +1** makes the original fast numerator assertion **FAIL**; restoring it makes the same checker pass (`method.log`).

   The fast method is sound: each denominator is a nonzero content times monic polynomial factors; their lcm clears all fractions while retaining contents, and each division has an asserted zero remainder. A zero resulting polynomial is an exact rational identity. Denominator clearing loses no admissible key.

5. **OK — E3, normalization, and domains; denominator wording needs correction.**

   Both E3 certificates pass at depths 9 and 10 (`E3_9.log`, `E3_10.log`): the direct DL2 and folded conjugate feeds cancel exactly. The diagonal normalization is anti-real, so Dc=-conj(Delta), as implemented. E3 proves blindness of the adjacent child; it does not place the initial diagonal coefficient on its stable polynomial.

   Depths 9,10 require d>=7. Adjacent d=2p-1 requires p>=4; diagonal d=2p likewise requires p>=4 here. Every feed denominator 2d-dp-1 is at least 3, and the other relevant denominators are nonzero. They are **not all odd**: odd dp gives an even feed denominator.

   The off-diagonal tops -16i(d-6)/Pi_4 and -61/Pi_5 are nonzero for d>=7. The stable diagonal tops, normalized by rho_pp(0), simplify independently to

   ```text
   depth 9:  (-32-32i)(p-3) / [(4p-9)(4p-7)(4p-5)(4p-3)]
   depth 10: (-61+61i) / [(4p-11)(4p-9)(4p-7)(4p-5)(4p-3)].
   ```

   Both are nonzero for p>=4 (`details.log`). Thus the new tower degrees are exactly floor(delta/2).

6. **OK — the level k=7 and the all-k weight-seven consequence.**

   WTS Lemma 3 is conservatively stated for k>=d+2, but its proof uses the fixed real-linear parent-to-kernel map at the same k. For kernel degree d=6, all degree-7 parents are **(7,0),(6,1),(5,2),(4,3)**. They are off-diagonal, so their proved tower law starts at k=7. Also j=i-1 is impossible when i+j=6, eliminating the extra folded diagonal-parent term. Since k=7>6, there is no initial work-item term. The same proof therefore covers k=7. In TOW4, a degree-6 diagonal parent forced a separate k=6 boundary check; that obstruction is absent here. This does not deny that degree-7 adjacent towers themselves require E3 when fed from degree 8; those identities were proved above.

   Independent `check_kernels.py` imports no submitted verifier. Exact Newton differences from k=8..13 reproduce **all 56 occupied real-pole QSI coordinates at every k=7..13: 392 equalities**, all four keys (`kernels.log`). Missing pole entries represent zero; the union of all poles is checked. Every +/-i part and constant part vanishes separately for every key at k=6..13. At k=6, **18/56** real coordinates differ from the quintic. The submitted `verify_W7.py` independently reruns successfully (`submitted_W7.log`). The number 56 already includes all four keys.

   Degree-7 towers through depth 10 exhaust their higher parts. The fixed real-linear kernel map gives degree-at-most-five residue polynomials for every k>=7. Their six +/-i zeros at k=8..13 force identical vanishing. The rational feed is proper (integrated towers are bounded at infinity and DL decays), so no polynomial part is missing. Applying the fixed integration map preserves polynomial degree. The six forward-difference kernels based at k=7 are saved in `W7_Newton_kernels.pkl`; define E_{6,m} as their fixed integral images. Consequently

   ```text
   W_7(k) = sum_{m=1}^6 binom(k-7,m-1) E_{6,m}, every integer k>=7.
   ```

   The k=6 discrepancy is proved here at the **kernel level**; no separate numerical inequality between integrated constants was evaluated. Finite kernel interpolation confirms the boundary but is not the all-k proof.

7. **Task E — OK for the finite non-adjacent reduction; GAP if promoted to a full all-depth equivalence; prefactor prose ERROR, readily repairable.**

   With P=floor(delta/2), U_dp=sum_s(-1)^s F(dp,s), and closing residual C_delta=Lambda_delta(P)-F(delta,P), finite telescoping gives

   ```text
   C_delta = (-1)^P [feed_delta(U)-a_delta].
   ```

   Both sides vanish exactly at delta=9,10 (`details.log`). Evaluation at x=-1 is legitimate coefficientwise for these finite Newton polynomials. It is not automatically defined for an arbitrary nonterminating formal x-series. A full all-depth tower reduction must additionally retain adjacency and E3 at every depth; identifying a separately prescribed zigzag-descending family also requires actual top matching, not just top equations. The d^5/Pi_5 perturbation in item 1 supplies a concrete counterexample to identification from closing/top equations alone. Those extra premises are separately established through depth 10, so the theorem stands.

   The submitted `taskE_defect.py` reports all exact for delta=2..10 (`submitted_defect.log`). Independently, taking D_delta to be the shared reduced a_delta denominator, the quotient in

   ```text
   U_delta-a_delta = prefactor_delta * Q_delta / D_delta
   ```

   is polynomial in (d,m) over Q(i), with degree_d Q_delta=3(floor(delta/2)-1), for every delta=2..10 (`details.log`). The prefactor is the stated even scalar or odd scalar times i(d-(delta+3)/2). Equal normalized denominators and the first three leading numerator coefficients at 9,10 also pass (`shapes.log`).

   `result.md` omits `/D_delta` while calling Q_delta polynomial; the later prose in the proof drops it again. That is false as written and must retain the denominator. Nonzero numerical scalar/sign factors alone have no divisibility content over Q(i); the odd linear factor is the substantive observation. These are finite-depth observations, not a universal theorem or an all-depth proof.

8. **Wording and reproducibility fixes.**

   - State the exact admissible domains and the separate initial diagonal exception. Replace “at k=2p E3 gives the stable extrapolation” with the child-blindness statement in item 3.
   - Replace “nonzero odd” denominators by “nonzero” denominators. Replace unrestricted “seed-blind” by “closing does not independently certify the seed,” or specify additive-constant perturbations.
   - Qualify claims that closing plus E3 alone is equivalent to the full folded tower law: the adjacent identities are additional requirements. Apply the qualifications in item 7 to the all-depth outlook.
   - Retain `/D_delta` in every prefactor formula and say “verified for delta=2..10.” A polynomial quotient per checked depth is not a general closed form for Q_delta.
   - All six fully displayed new F formulas and the full printed N_5 are correct by mechanical symbolic comparison (`details.log`). Both a_9/a_10 pickle formulas also pass. Point to those pickles for elided numerators.
   - Correct the evidence count: **56 total QSI coordinates across four keys**, not four times 56. Use degree **at most five** for the integrated law unless separately asserting its leading coefficient is nonzero.
   - `verify_tow5.py` checks child return codes, but its own final summary does not exit nonzero when a reported lemma fails. `forced910.py` and `taskE_defect.py` likewise should assert their final success flags. The official runner also hard-codes `nice`, and the seed/cache writers target their own folder. Provide an output-directory option and direct subprocess execution. These are reproducibility repairs; the actual residuals were inspected here.
   - Update stale library comments calling R_3/R_4 fitted or pending and the R5_raw provenance reference in `tow5lib.py`.

**Fresh confirmation on the supplied full-tower data.** Both scripts pass: `check_full910.py` gives **990 stable cells (132 diagonal), zero bad, plus 11 initial diagonal cells, zero bad**; `check_full910_ext.py` gives **16,225 stable cells (1,122 diagonal), zero bad, plus 66 initial diagonal cells, zero bad**. See `check_full910.log`, `check_full910_ext.log`. These are checks of the supplied sparse data, not a regeneration of all those cells.

**Execution notes and reproduction.** Local copies preserve the submitted inputs; only the local forced cache was removed to ensure a rebuild. `run_certificate.py` substitutes exact Gaussian-rational field cancellation for generic cancellation; for `verify_L.py` it also pre-cancels the same cumulant expressions. The original equations, variable domains, and zero tests remain unchanged. The raw L and forced runs were interrupted after the partial results explicitly described above because of expression growth. The first exploratory details run contained my incorrect assertion that a constant seed perturbation must be detected; it failed (`details_first.log`). I corrected that assertion to d^6/Pi_5, then reran the entire details check with exit 0. This was a review-code error, not a failed submitted certificate. The separate `check_seed.py` now checks both the invisible constant and detected degree-six perturbations.

Use the specified venv Python with `-B`, cwd this folder. After setup, the following reproduce the checks; `setup_review.py` removes the local forced cache, so rerun forced before its consumers:

```sh
PY='<project>/.venv/bin/python'
"$PY" -B setup_review.py
"$PY" -B prove_S5.py
"$PY" -B check_r5.py
"$PY" -B run_certificate.py verify_L.py
"$PY" -B run_certificate.py forced910.py 10
"$PY" -B run_certificate.py check_initial.py
"$PY" -B run_adjacent.py
"$PY" -B run_certificate.py check_method.py
"$PY" -B run_certificate.py prove_E3b910.py 9
"$PY" -B run_certificate.py prove_E3b910.py 10
"$PY" -B check_kernels.py
"$PY" -B verify_W7.py
"$PY" -B run_full.py
"$PY" -B run_certificate.py check_details.py
"$PY" -B check_shapes.py
"$PY" -B run_certificate.py check_seed.py
"$PY" -B run_certificate.py ../TOW5_tower_depth10/taskE_defect.py
```

**Limits.** I did not run the monolithic `verify_tow5.py`, the optional mpmath tau cross-check, the complete 425-per-depth stored-rho loop, all inherited base certificates, or a fresh reduction generating the WTS kernels. I did independently generate the 132 initial and 20 propagated tower coefficients specified above. The corrected depth<=8 theorem, E0/T1 reduction, and written WTS fixed-map argument are inherited premises. No source, data, dependency, memory, or paper outside this folder was modified. Depth>=11, an all-depth closure proof, and weight>=8 remain unproved by this review.
