# TOWN1 — exact test of the tower law T at LARGE degree n (n ≤ 40, data to n ≤ 46), small depth (δ ≤ 8)

## VERDICT: **T holds in every tested cell, with zero exceptions. No counterexample (Task D) exists.**

The tower law **T** (WTS-theory §2.3 / TOW1) — for a key `(p,q)` of degree `n = p+q`, the depth-δ tower
coefficient `u^{(k)}_{pq}(δ)` (the coefficient of `(v−i)^{−(2n−2−δ)}` in the higher part of `W^{(k)}_{pq}` at
`v = +i`) is, for `k ≥ n` (diagonal key `(n/2,n/2)` of even `n`: `k ≥ n+1`), a polynomial in `k` of degree
**exactly** `floor(δ/2)` — is confirmed **exactly** in `Q(√2, i)`:

| depth cap D | degrees n | K (max level) | genuine cells | exceptions | min check-points/cell |
|---|---|---|---|---|---|
| 6 | 2..40 | 52 | **3047** | 0 | 8 |
| 8 | 2..40 | 52 | **3905** | 0 | 7 |
| 8 | 2..45 | 52 | **4904** | 0 | 3 |
| 8 | 2..46 | 52 | 5120 | 0 | 1 |

Every cell is a *genuine* test (`#levels ≥ floor(δ/2)+2`, ≥1 exact check point beyond the fit). Test: the
`(floor(δ/2)+1)`-th finite difference of the exact level sequence vanishes identically (exact QSI zeros)
**and** the `floor(δ/2)`-th is nonzero (degree exactly `floor(δ/2)`, leading coefficient nonzero).
**Diagonal shift confirmed exactly**: diagonal `(p,p)` value at `k = n = 2p` is *off* the polynomial in all 110
diagonal cells probed; off-diagonal `k = n` value is *on* it in all 3736 off-diagonal cells probed.

Evidence for T (hence, via WTS Lemma 3, for the kernel law K(d)) is pushed from `n ≤ 11, k ≤ 14` (64 cells,
WTS-theory) to **`n ≤ 40 (46), k ≤ 52`, depths `δ ≤ 8`** (thousands of cells).

## Method and why it is exact (the truncated engine, Task A)

By WTS Lemma 4 / TOW1 §2, `Λ` is **depth non-decreasing**, so the depth-≤D part of the tower is determined by
the depth-≤D part of the initial data `ρ_m` and of `Λ`. Both are computed from Laurent/Taylor series at
`v = i` (series in `s = v−i`, truncated at `s^D`), never building the full integrand PF.

Initial data (proved formula; the AB-part has no pole at `v=i`), key `(p,q)`, degree `d`:
`Φ(v) = −(√2/2)·8^d·v^{2d−1}·S(t(v))/[(v²−a²)(v²−b²)(v²+1)^{2d−2}]`, `t(v)=(v²+2v−1)/(4v)`, and
`ρ_{pq}(δ) = [s^δ] Common(s)·[S_p(t(i+s)) + [p>q]·S_q(1−t(i+s))]` with
`Common(s) = −(√2/2)·8^d·(i+s)^{2d−1}·(2i+s)^{2−2d}/G(i+s)`, `G = v⁴−6v²+1`. The operator is TOW1's proved
Prop. T1 (`Lambda_ij`), truncated to `δ ≤ D`. Everything is exact `QSI = Q(√2,i)` (4 `Fraction`s), so all
checks are exact equalities (residual identically 0, not a numerical approximation).

**Task A validation (identical coefficients to independent data):**

| check | coefficients | mismatches (residual) |
|---|---|---|
| ρ vs `towlib.rho`, d=2..14, δ≤min(6,2d−4) | 408 | 0 (exact) |
| tower vs `towlib.u`, k=2..14, all keys, δ≤min(6,2n−4) | 2083 | 0 (exact) |
| tower vs `towlib.u`, k=2..14, all keys, δ≤min(8,2n−4) | 2543 | 0 (exact) |
| ρ vs full integrands `binteg_15..18` (independent build), δ≤min(8,2d−4) | 324 | 0 (exact) |

## Task C — consequences (EXACT for the stated k-ranges)

Kernels computed from the exact tower values via WTS Lemma 2 (`kernel_from_tower.py`); construction validated
against every recorded kernel (WT3 k=3..15, WTS `kern_K` k=4..18, WT4 k=4..12: all identical).

* **Degree-2 kernels (weight 3) R_20(k;v), R_11(k;v) — EXACT for k = 3..52**, matching the WT3 closed forms
  (0 mismatches), extending WT3's exact `k ≤ 15` to **`k ≤ 52`**:
  `R_20 = (1−√2)/(4v) + 1/(4(v−a)) + (2k−7)/(4(v+a)) − (k−2)/(2(v+b))`,
  `R_11 = (2k−5+√2)/(4v) − (2k−5)/(2(v+a))`. Hence (WT3's proved reduction) `W3(k) = W3′ − (k−3)S` exact for
  `3 ≤ k ≤ 52`.
* **Degree-3 kernel law (weight 4) — EXACT for k = 4..52**: every real-pole residue of `K_30`, `K_21` is a
  polynomial in `k` of degree ≤ 2 (3rd finite differences vanish exactly), each key degree exactly 2, no `±i`
  pole survives (explicit polynomials in `results/kernel_laws.json`). Since WT4's only `k`-dependent input is
  this kernel law, **`W_4(k) = W4″ + (k−4)D1 + binom(k−4,2)D2` is now exact for `4 ≤ k ≤ 52`** (was `k ≤ 14`).

## Deliverables

- `result.md` (this result, self-written), `progress.log` (one line/step).
- `trunc_engine.py` (engine), `driver.py` (harvest/test), `kernel_from_tower.py` (Lemma-2 kernels).
- `data/towers_D8_K38.pkl`, `data/towers_D8_K48.pkl`, `data/towers_D8_K52.pkl` (harvested exact `u^{(k)}_{pq}(δ)`).
- `results/cells_D6_N40.json`, `results/cells_D8_N40.json`, `results/cells_D8_N45.json` (per-cell verdicts),
  `results/kernel_laws.json` (exact kernel polynomials).
- `verify_town1.py`: standalone re-verification, **5/5 PASS in ~3–20 s** (< 15 min): A1 ρ==towlib.rho (d=2..12,
  303 coeffs); A2 tower==towlib.u (k=2..12, 1323 coeffs); B1 tower law all genuine cells n=2..18 at D=6 (659
  cells, 0 fails); B2 diagonal shift (41/41 off); C1 R_20,R_11==closed forms at k=16..22. Re-confirmed PASS
  after the interruption.

## Evidence / step files

All arithmetic is exact (`QSI` over `Q(√2,i)`); every reported check has residual **exactly 0** (integer
equality, not a float tolerance — hence no "digit count" applies; the evidence unit is exact coefficients
compared, tabulated above, all with 0 mismatches). The computation is recorded in the harness-saved
`step_01.py … step_56.py` in this folder; the computational chain is: `step_01/02` (ρ engine + validation),
`step_03–05` (tower build + validation vs `towlib.u`, D=6,8), `step_06/07` (harvest D=8 K=38→48), `step_08/09`
(tower-law test N=40 + diagonal shift), `step_12–16` (Task C: degree-2/3 kernels vs WT3/WTS/WT4 records),
`step_17–19` (extend to K=52, retest, re-export cells), `step_20` (Task C to k=52), `step_22` (independent ρ
check vs `binteg_15..18`), `step_23` (final verify 5/5 PASS). Reusable modules: `trunc_engine.py`,
`driver.py`, `kernel_from_tower.py`, `verify_town1.py`. `fastred.py`/`genintegrand.py` are imported read-only
from `WTS_weight_k_structure/theory/`; nothing was written into other job folders.