import os
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)

md = r'''# TOWN1 — exact test of the tower law T at LARGE degree n (n ≤ 40), small depth (δ ≤ 8)

Folder: `research/08_theory/TOWN1_tower_large_n/`.  Python: project venv (sympy/mpmath/python-flint available;
only `sympy` and the stdlib `fractions` are used, via the copied `fastred.py` QSI arithmetic).
All scripts run one at a time under `nice -n 10`, each well under 30 min.  Nothing written outside this folder.

## VERDICT

**T holds in every tested cell, with zero exceptions.**  The tower law T of WTS-theory §2.3 / TOW1 —
for a key `(p,q)` of degree `n = p+q`, the depth-δ tower coefficient `u^{(k)}_{pq}(δ)` (coefficient of
`(v−i)^{−(2n−2−δ)}` in the higher part of `W^{(k)}_{pq}` at `v = +i`) is, for `k ≥ n` (for the diagonal key
`(n/2,n/2)` of even `n`: `k ≥ n+1`), a polynomial in `k` of degree **exactly** `floor(δ/2)` — is confirmed
**exactly** (in `Q(√2, i)`) for:

| depth cap D | degrees n | cells tested (genuine) | exceptions |
|---|---|---|---|
| 6 | 2..30 | 1752 | 0 |
| 8 | 2..30 | 2240 | 0 |
| 6 | 2..40 | 3047 | 0 |
| 8 | 2..40 | 3905 | 0 |

Every one of these cells is a *genuine* test (number of levels `≥ floor(δ/2)+2`, i.e. at least one exact
check point beyond the `floor(δ/2)+1` points that determine the polynomial; minimum over all cells is
3–4 check points at `D=8, N=40`).  For each cell the test is: the `(floor(δ/2)+1)`-th finite difference of
the level sequence vanishes identically (exact `QSI` zeros) **and** the `floor(δ/2)`-th finite difference is
nonzero (degree exactly `floor(δ/2)`, leading coefficient nonzero).

**Diagonal shift confirmed exactly.**  For the diagonal keys `(p,p)` the value at the initial level `k = n = 2p`
is **off** the stable polynomial in all 110 diagonal cells probed (the law starts at `k = 2p+1`), while for
off-diagonal keys the `k = n` value is **on** the polynomial in all 3736 off-diagonal cells probed (`k ≥ n`).

This pushes the exact evidence for T (hence, via TOW1/WTS Lemma 3, for the kernel law K(d)) from the
previous `n ≤ 11, k ≤ 14` (64 cells, WTS-theory) up to **`n ≤ 40, k ≤ 48`, depths `δ ≤ 8`** (thousands of cells).

## 0. Method — the truncated engine (Task A)

The point (WTS Lemma 4 / TOW1 §2): the recursion operator `Λ` is **depth non-decreasing**, so the depth-≤ D
part of the whole tower is determined by the depth-≤ D part of the initial data `ρ_m` and of `Λ`.
We compute both directly from Laurent/Taylor series at `v = i` (series in `s = v − i`, truncated at `s^D`),
**never** building the full integrand PF.

*Initial data (proved formula, TOW1 §3; the AB-part has no pole at `v=i`).*  For key `(p,q)`, degree `d`:
```
Φ(v) = −(√2/2)·8^d · v^{2d−1} · S(t(v)) / [(v²−a²)(v²−b²)(v²+1)^{2d−2}],    t(v) = (v²+2v−1)/(4v)
ρ_{pq}(δ) = [s^δ]  Common(s)·[ S_p(t(i+s)) + [p>q]·S_q(1−t(i+s)) ]
Common(s) = −(√2/2)·8^d·(i+s)^{2d−1}·(2i+s)^{2−2d} / G(i+s),   G(v) = v⁴−6v²+1 = (v²−a²)(v²−b²)
S_p(t(i+s)) = (1/(p!q!)) Σ_j C_j Σ_s S(p,s)S(q,j−s) M[(j,s)],   C_j = binom(1/2,j)·j!·2^d
M[(j,s)] = T^{2(j−s)} U^{2s} W^{d−j},   T = t(i+s), U = 1−T, W = q₂(t(i+s)) = T²+U²  (series mod s^{D+1})
```
(the `4^{−D} v^{2d−1−D}` and the `(4v)^D` from `S̃(v) = (4v)^D S(t(v))` cancel exactly, for any `deg S`).
The series `M[(j,s)]` are shared by all keys of the level, so one level costs `O(d²D²)` setup + `O(d³D)` accumulation.

*Operator (TOW1 Prop. T1, proved).*  `u^{(k)}_{pq}(δ) = Λ_{ij}` of the degree-`(n+1)` data at the same level `k`
(the two `τ1/τ2` parent sums plus the diagonal fold term), truncated to `δ ≤ D`.  τ-coefficients and `i`-powers
are memoised; the whole engine is exact `QSI = Q(√2,i)` (4 `Fraction`s).

**Task A validation (identical coefficients to the recorded towers `towlib.u`, `towlib.rho`, k ≤ 14):**

| check | coefficients compared | mismatches |
|---|---|---|
| ρ vs `towlib.rho`, d = 2..14, δ ≤ min(6,2d−4) | 408 | 0 |
| tower vs `towlib.u`, k = 2..14, all keys, δ ≤ min(6,2n−4) | 2083 | 0 |
| tower vs `towlib.u`, k = 2..14, all keys, δ ≤ min(8,2n−4) | 2543 | 0 |

## 1. Task B — the cell table

Harvest `u^{(k)}_{pq}(δ)` for **every** key and every level `k = 2..K`, depths `δ ≤ min(D,2n−4)`
(`data/towers_D8_K38.pkl`, `data/towers_D8_K48.pkl`), then the finite-difference test per cell
(`driver.test_tower_law`).  A cell is *genuine* iff `#levels ≥ floor(δ/2)+2`.

Per-degree cell counts at `D=8, N=40` (all pass): n=2:2, n=3:6, n=4:15, n=5:21, then for n ≥ 6 with `δ ≤ 8`
capped by `2n−4`: n even/odd give `9·(keys of degree n)`; every degree n=2..40 has full depth coverage `δ = 0..8`.
At `D=6, N=40`: 3047 cells; at `D=8, N=40`: 3905 cells.  Full machine-readable table: `results/cells_*.json`.

Timings (machine heavily loaded): engine build/validate ≤ 14 ≤ 4 s; harvest D=8 K=2..38: 67 s;
extend K=39..48: 71 s; per level at k≈48: ~9 s (ρ ~2.5 s + tower ~6 s).

## 2. Task C — consequences (EXACT for k ≤ K)

(see below; filled in by step_C.)

## 3. Files

| file | purpose |
|---|---|
| `trunc_engine.py` | the truncated engine: series arithmetic, `level_series`, `rho_key`, `rho_level`, `Lambda_ij` (T1), `build_tower_level` |
| `driver.py` | `harvest` (incremental, resumable), `fd_test` (finite-difference polynomial test), `test_tower_law` |
| `data/towers_D8_K38.pkl`, `data/towers_D8_K48.pkl` | harvested exact tower coefficients `u^{(k)}_{pq}(δ)` |
| `results/cells_*.json` | the per-cell verdicts |
| `verify_town1.py` | reruns Task A validation (k ≤ 12) and a sample of Task B cells, PASS/FAIL, < 15 min |
| `progress.log` | one line per finished step |

Reuse note: `fastred.py`/`genintegrand.py` are **imported** (read-only) from `WTS_weight_k_structure/theory/`;
nothing is written into other job folders.
'''
with open('result.md','w') as f:
    f.write(md)
print("wrote result.md", len(md), "bytes")
with open('progress.log','a') as f:
    f.write("step_10 wrote result.md (Tasks A,B)\n")
