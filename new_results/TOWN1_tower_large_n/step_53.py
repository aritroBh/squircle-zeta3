import os
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
# clean temp verify pickle (verify_town1.py regenerates it)
tmp = 'data/_verify_towers_D6_K22.pkl'
if os.path.exists(tmp): os.remove(tmp)

md = r'''# TOWN1 — exact test of the tower law T at LARGE degree n (n ≤ 40, data to n ≤ 46), small depth (δ ≤ 8)

Folder: `research/08_theory/TOWN1_tower_large_n/`.  Python: project venv (only `sympy` and stdlib `fractions`
are used, via the read-only imported `fastred.py` QSI arithmetic from `WTS_weight_k_structure/theory/`).
All scripts run one at a time under `nice -n 10`, each well under 30 min.  Nothing written outside this folder.

## VERDICT

**T holds in every tested cell, with zero exceptions.**  The tower law T of WTS-theory §2.3 / TOW1 —
for a key `(p,q)` of degree `n = p+q`, the depth-δ tower coefficient `u^{(k)}_{pq}(δ)` (the coefficient of
`(v−i)^{−(2n−2−δ)}` in the higher part of `W^{(k)}_{pq}` at `v = +i`) is, for `k ≥ n` (for the diagonal key
`(n/2,n/2)` of even `n`: `k ≥ n+1`), a polynomial in `k` of degree **exactly** `floor(δ/2)` — is confirmed
**exactly** (in `Q(√2, i)`) for:

| depth cap D | degrees n | K (max level) | genuine cells | exceptions | min check-points/cell |
|---|---|---|---|---|---|
| 6 | 2..40 | 52 | **3047** | 0 | 8 |
| 8 | 2..40 | 52 | **3905** | 0 | 7 |
| 8 | 2..45 | 52 | **4904** | 0 | 3 |
| 8 | 2..46 | 52 | 5120 | 0 | 1 |

Every cell is a *genuine* test: the number of levels is `≥ floor(δ/2)+2`, i.e. there is at least one exact
check point beyond the `floor(δ/2)+1` points that determine the polynomial.  The test is: the
`(floor(δ/2)+1)`-th finite difference of the exact level sequence vanishes identically (exact `QSI` zeros)
**and** the `floor(δ/2)`-th finite difference is nonzero (degree exactly `floor(δ/2)`, leading coefficient
nonzero).  At the headline `D=8, N=40` every cell has ≥ 7 independent check points; the data further covers
all degrees to `n ≤ 46` at `D=8`.

**Diagonal shift confirmed exactly.**  For the diagonal keys `(p,p)` the value at the initial level `k = n = 2p`
is **off** the stable polynomial in all 110 diagonal cells probed (the law starts at `k = 2p+1`), while for
off-diagonal keys the `k = n` value is **on** the polynomial in all 3736 off-diagonal cells probed.

This pushes the exact evidence for T (hence, via WTS Lemma 3, for the kernel law K(d)) from the previous
`n ≤ 11, k ≤ 14` (64 cells, WTS-theory) up to **`n ≤ 40 (46), k ≤ 52`, depths `δ ≤ 8`** (thousands of cells).
No cell of Task D (a failing cell to characterise) exists.

## 0. Task A — the truncated engine and its validation

The point (WTS Lemma 4 / TOW1 §2): the recursion operator `Λ` is **depth non-decreasing**, so the depth-≤ D
part of the whole tower is determined by the depth-≤ D part of the initial data `ρ_m` and of `Λ`.  Both are
computed directly from Laurent/Taylor series at `v = i` (series in `s = v − i`, truncated at `s^D`), never
building the full integrand PF.

*Initial data (proved formula, TOW1 §3; the AB-part has no pole at `v = i`).*  For key `(p,q)`, degree `d`:
```
Φ(v) = −(√2/2)·8^d · v^{2d−1} · S(t(v)) / [(v²−a²)(v²−b²)(v²+1)^{2d−2}],   t(v) = (v²+2v−1)/(4v)
ρ_{pq}(δ) = [s^δ]  Common(s)·[ S_p(t(i+s)) + [p>q]·S_q(1−t(i+s)) ]
Common(s) = −(√2/2)·8^d·(i+s)^{2d−1}·(2i+s)^{2−2d} / G(i+s),   G(v) = v⁴−6v²+1 = (v²−a²)(v²−b²)
S_p(t(i+s)) = (1/(p!q!)) Σ_j C_j Σ_s S(p,s)S(q,j−s) M[(j,s)],   C_j = binom(1/2,j)·j!·2^d
M[(j,s)] = T^{2(j−s)} U^{2s} W^{d−j},   T = t(i+s), U = 1−T, W = q₂(t(i+s)) = T²+U²   (series mod s^{D+1})
```
(the `4^{−D} v^{2d−1−D}` and the `(4v)^D` of `S̃(v) = (4v)^D S(t(v))` cancel exactly, for any `deg S`; the
series `M[(j,s)]` are shared by all keys of the level, so one level costs `O(d²D²)` setup + `O(d³D)`).
The operator is TOW1's Prop. T1 (`Lambda_ij`), truncated to `δ ≤ D`; τ-coefficients and `i`-powers memoised;
everything exact `QSI = Q(√2,i)` (4 `Fraction`s).  Implementation: `trunc_engine.py`.

**Validation against the recorded towers (`towlib.rho`, `towlib.u`), identical coefficients:**

| check | coefficients compared | mismatches |
|---|---|---|
| ρ vs `towlib.rho`, d = 2..14, δ ≤ min(6,2d−4) | 408 | 0 |
| tower vs `towlib.u`, k = 2..14, all keys, δ ≤ min(6,2n−4) | 2083 | 0 |
| tower vs `towlib.u`, k = 2..14, all keys, δ ≤ min(8,2n−4) | 2543 | 0 |

## 1. Task B — the cell table

Harvest `u^{(k)}_{pq}(δ)` for **every** key and every level `k = 2..52`, depths `δ ≤ min(8,2n−4)`
(`data/towers_D8_K52.pkl`, also `K38`, `K48`), then the finite-difference test per cell (`driver.test_tower_law`).
Per-degree cell counts at `D=8` (all pass): n=2:2, n=3:6, n=4:15, n=5:21, and for n ≥ 6 every degree has full
depth coverage `δ = 0..8` (capped by `2n−4`), i.e. `9 × (#keys of degree n)` cells.  Machine-readable tables:
`results/cells_D6_N40.json`, `results/cells_D8_N40.json`, `results/cells_D8_N45.json` (each with per-cell
`npts`, `nchecks`, `pass`).

Timings (machine heavily loaded): engine build+validate (k ≤ 14) ≈ 4 s; harvest D=8 K=2..38: 67 s;
K=39..48: 71 s; K=49..52: 55 s (≈ 13 s/level at k ≈ 50).

## 2. Task C — consequences (EXACT for the stated k-ranges)

Kernels are computed from the **exact tower values** via WTS Lemma 2 (`kernel_from_tower.py`);
the construction was validated against every recorded kernel (WT3 `wt3_kernels.json` k=3..15, WTS `kern_K.pkl`
k=4..18, WT4 `kernels_deg3.json` k=4..12): all identical.

* **Degree-2 kernels (weight 3), R_20(k;v), R_11(k;v) — EXACT for k = 3..52.**  Computed from the degree-3
  tower (`W_20 = −3 Gm_30 DL1 − Gm_21 DL2`, `W_11 = −2 Gm_21 DL1`, simple part), they match the WT3 closed
  forms **exactly** for every `k = 3..52` (0 mismatches):
  ```
  R_20(k;v) = (1−√2)/(4v) + 1/(4(v−a)) + (2k−7)/(4(v+a)) − (k−2)/(2(v+b))
  R_11(k;v) = (2k−5+√2)/(4v) − (2k−5)/(2(v+a)),        a = 1+√2, b = √2−1
  ```
  This extends the exact kernel verification from WT3's `k ≤ 15` to **`k ≤ 52`** (and, via WT3's proved
  reduction theorem, `W3(k) = W3′ − (k−3)S`, `S = ζ(3)/8 + π²log2/24 + log³2/12`, is exact for `3 ≤ k ≤ 52`).

* **Degree-3 kernel law (weight 4) — EXACT for k = 4..52.**  Computed from the degree-4 tower
  (`W_30 = −4 Gm_40 DL1 − Gm_31 DL2`, `W_21 = −3 Gm_31 DL1 − 2 Gm_22 DL2 − 2 F(Gm_22 DL1)`; note the diagonal
  fold).  Every real-pole residue of `K_30`, `K_21` is a polynomial in `k` of degree ≤ 2 across `k = 4..52`
  (3rd finite differences vanish exactly), each key has degree exactly 2, and no `±i` pole survives:
  ```
  K_30(k;v): (v+a): (1/6)k² − (2/3)k + 5/24 ;  v: −(1/12)k + 3/8 − √2/8 ;
             (v−a): −(1/6)k² + k − 35/24 ;   (v+b): −(1/6)k + 1/2
  K_21(k;v): (v+a): −(1/2)k² + (5/2)k − 23/8 ;  v: (1/4)k² − (5/4)k + 23/16 + √2/8 ;
             (v−a): (1/4)k² − (3/2)k + 37/16 ; (v+b): −(1/4)k² + (3/2)k − 37/16
  ```
  (coefficient of `1/(v−pole)` = the displayed polynomial in `k`; exact for `4 ≤ k ≤ 52`, see
  `results/kernel_laws.json`).  Since the only `k`-dependent input to WT4's formula
  `W_4(k) = W4″ + (k−4)·D1 + binom(k−4,2)·D2` (with `W4″`, `D1`, `D2` proved constants) is the degree-3
  kernel polynomial law, **this formula is now exact for `4 ≤ k ≤ 52`** (was `k ≤ 14`).

## 3. Files

| file | purpose |
|---|---|
| `trunc_engine.py` | truncated engine: QSI series arithmetic, `level_series`, `rho_key`, `rho_level`, `Lambda_ij` (T1), `build_tower_level` |
| `driver.py` | `harvest` (incremental, resumable), `fd_test`, `test_tower_law` |
| `kernel_from_tower.py` | Lemma-2 tower→kernel: `gm_from_tower`, `kernel_20_21_level`, `kernel_deg3_level`, closed forms `R20_closed`, `R11_closed` |
| `data/towers_D8_K38.pkl`, `data/towers_D8_K48.pkl`, `data/towers_D8_K52.pkl` | harvested exact tower coefficients `u^{(k)}_{pq}(δ)` |
| `results/cells_D6_N40.json`, `results/cells_D8_N40.json`, `results/cells_D8_N45.json` | per-cell verdicts |
| `results/kernel_laws.json` | the exact kernel polynomials and verification ranges |
| `verify_town1.py` | standalone re-verification (below), **5/5 PASS in 20 s** |
| `progress.log` | one line per finished step |

`verify_town1.py` reruns: A1 ρ == `towlib.rho` (d=2..12, 303 coeffs); A2 tower == `towlib.u` (k=2..12, all keys,
1323 coeffs); B1 tower law T all genuine cells n=2..18 at D=6 (harvest K=22, 659 cells, 0 fails); B2 diagonal
shift sample (41/41 off-polynomial); C1 degree-2 kernels == closed forms R_20,R_11 at k=16..22.  Prints
PASS/FAIL per check; final `ALL PASS`.

Reuse note: `fastred.py`/`genintegrand.py` are **imported** (read-only) from `WTS_weight_k_structure/theory/`;
the recorded towers/kernels are read from `WTS`, `WT3`, `WT4`; nothing is written into other job folders.
'''
with open('result.md','w') as f:
    f.write(md)
print("wrote result.md", len(md), "bytes")
with open('progress.log','a') as f:
    f.write("step_21 result.md finalized (Tasks A,B,C, verify)\n")
