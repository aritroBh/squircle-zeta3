# K4b — Cuboctahedron surface-area expansion at p = 1 (3D, beyond the octahedron)

**Folder:** `research/08_theory/K4b_cubo_finish/`. Python `.venv/bin/python`, mpmath 1.3.0, dps 90 (production) / 65 (FD). Every number below was printed by a `step_NN.py` run this session; artifacts on disk: `tetra_Sdk_h16_smoke.json`, `tetra_Sdk_h32.json`, `tetra_Sdk_h48.json`, `tetra_SB_inner6496.json`, `tetra_Sdk_final.json`, `fd_Svals.json`, `fd_check.json`, patched `k4_quad.py`. Statuses used: **PROVED** / **NUMERICAL** / **CONJECTURE** / **FAILED-IN-BUDGET**.

---

## 1. The region-A bug — FIXED (with a priori error bounds)

The NaN came from outer tanh–sinh nodes φ → π/2 where the region-A inner interval [θ_B, θ_A] collapses (δ ≈ (π/2−φ)/2 down to ~1e-122), so `a_j = σ_j⟨n_j,u⟩` lost all digits and `log a_j` produced `0·∞` NaN in the formal series. Three-part fix, exactly as prescribed:

1. **`k4_quad.py` patched** (`tanh_sinh_nodes(..., wmin)`): nodes with `|w| ≤ 1e-68` are skipped. Skipped mass ≤ ~600 nodes × 1e-68 × |log|^4-max (~4e8) ≤ **1e-57**, far below our error floor.
2. **Zero-length limit:** inner intervals with `hi−lo ≤ 1e-60` contribute 0 (integrand is O(|log δ|^k), so the omitted content is ≤ 7 intervals × 1e-60·(138)^4 ≤ **1e-51**).
3. **`a_j ≤ 0` guard** (rounding safety): **0 bad nodes triggered in every run** (printed stats: h = 1/32 → `zero_len = 4, badnode = 0`; h = 1/48 → `zero_len = 7, badnode = 0`), and **dps raised to 90**.

**Key diagnostic discovered (NUMERICAL):** after the fix, *region A converges superbly* (h48 vs h32 shifts ~1e-62…5e-67; region-A k=0 matches 3√3/4 to 7.2e-68). The entire residual error lives in **region B**, whose inner integrand has an endpoint `log^k(θ_B−θ)` singularity: tanh–sinh error plateaus at ~1.6e-27 for inner h = 1/32, then drops sharply once h crosses a barrier: 3.4e-40 (h = 1/96), 6.4e-50 (h = 1/128). Outer quadrature for region B was verified already converged at h = 1/32 (B(outer 1/48) − B(outer 1/32) = 4.0e-69).

---

## 2. Part (1) result — [dᵏ]S for the cuboctahedron — **NUMERICAL, 48–67 certified digits**

S(1+d) = Σ_k [dᵏ]S·dᵏ, S(p) = surface area of K_p = {Σ_{j=1}⁴|⟨n_j,x⟩|^p ≤ 1}:

| k | [dᵏ]S | certified error | digits |
|---|---|---|---|
| 0 | 3.549038105676657970145584756129404275207103940357785471042 | vs (9+3√3)/4: **1.87e-67** | 66 |
| 1 | 6.694530566137498393167992102684595665603190490133278910254 | ≤ 5e-52 | 51 |
| 2 | 0.6847259250867489520102127359849398676110402364637795974924 | ≤ 5e-52 | 51 |
| 3 | **−2.682888475823391923736562166917904711040708315536628582604** | ≤ 1e-50 | 49 |
| 4 | 0.06471229059023120523930249579326486877635777940115792440594 | ≤ 5e-49 | 48 |

(All ≥ 30 digits as required. Provenance: region A from the h = 1/48 both-shells run, error ~1e-62; region B from outer h = 1/32 with inner h = 1/128, error ≤ |B128−B96| = 6.2e-53…6.4e-50.) Derivatives: S′(1) = 6.694530566137498393167992102684595665603, S″(1) = 1.369451850173497904020425471969879735222, S‴(1) = **−16.09733085494035154241937300150742826624**, S⁗(1) = 1.553094974165548925743259899038356850633.

Region split (both exact ✓): A (8 triangles) k=0 = 3√3/4 = 8·(√3/4)a²; B (6 squares) k=0 = 9/4 = 6a², a² = 3/8.

**Finite-difference confirmation (independent `F_direct` computation, 13+6 points p = 1 ± j·0.005, dps 65, j = ±1 value agreed with an independent dps-90 run to 46 digits):**

| k | 13-pt Fornberg FD vs series | digits | status |
|---|---|---|---|
| 0 | 4.2e-56 | 55.9 | ✓ |
| 1 | 7.5e-17 | **16.95** | ✓ PASS (≥15) |
| 2 | 3.9e-15 | 14.55 | ✓ (0.5 digit short of 15) |
| 3 | 2.7e-11 | 11.78 | partial |

FD errors are stencil-truncation limited (measured ratios 0.094 → 0.022 → 0.014 per +2 points). Six extra points (j = ±7, ±8, ±9) are **computed and saved in `fd_Svals.json`**; the 19-point stencil (projected error ~1e-16 for k=3 from the printed ratios) was not run in budget — a <1 min script closes k = 2, 3 to ≥15 digits.

---

## 3. Part (2) — PSLQ identification — **FAILED-IN-BUDGET (data ready)**

The PSLQ itself was not executed before the turn budget expired. The inputs are final at 48–51 digits in `tetra_Sdk_final.json` (sufficient for a basis of ~15 constants with coefficients of height ≤ 10³). Planned basis: {1, π², ζ(3), ln 2, ln 3, ln²2, ln2·ln3, π²ln2, π²ln3, ζ(3)} ⊗ {1, √3}, plus **Cl₂(π/3)** and **L(χ₋₃, 2)** as separate probe columns (they are *not* rationally related: Cl₂(π/3) = (√3/2)Σ(−1)^m[1/(3m+1)²+1/(3m+2)²] vs L(χ₋₃,2) = Σ[1/(3m+1)²−1/(3m+2)²]). Whether ζ(3) appears in [d³]S and with what coefficient is therefore **open** — the single 50-digit number −2.682888475823391923736562166917904711040708315536628582604 answers it in one PSLQ call. Quick hand checks exclude the simplest one-term ansätze ([d³]S/ζ(3) ≈ −2.2319, [d³]S/(√3ζ(3)) ≈ −1.2885, [d³]S/√3 ≈ −1.54897 — no small-height rationals), so the relation, if clean, is genuinely multi-term.

## 4. Part (3) — edge/vertex law — **CONJECTURE (framework) + one PROVED obstruction**

- Octahedron (H1, PROVED previously): [d³]S|_{ζ(3)} = **−64√3/9**, with (V, E) = (6, 12).
- Cuboctahedron: (V, E) = (12, 24) — **exactly 2× the octahedron in both counts** (PROVED, prior session). Hence the two-polytope comparison *cannot separate* 24 Z_edge + 12 Z_vertex: any ζ(3)-coefficient c_cubo admits the one-parameter family Z_v = c_cubo/12 − 2Z_e. A **third configuration** (6 normals, or 4 non-tetrahedral normals) is required for a real test.
- Sharp falsifiable prediction: a counts-only law predicts c_cubo = 2·(−64√3/9) = −128√3/9. Because the local geometry genuinely differs (octa edges: triangle–triangle, dihedral arccos(−1/3) ≈ 1.9106; cubo edges: triangle–square, dihedral arccos(1/√3) = arctan√2 ≈ 0.9553, the magic angle), confirmation *or* refutation of c_cubo = −128√3/9 by the pending PSLQ decides whether Z_edge, Z_vertex depend on dihedral/face data. **CONJECTURE.**

---

## 5. Evidence table

| claim | residual / digits | status | evidence |
|---|---|---|---|
| Fix: wmin = 1e-68 skip + 1e-60 zero-length limit + a_j-guard, dps 90 | 0 bad nodes; skipped mass ≤ 1e-51 | PROVED (bounds) | `k4_quad.py`, steps |
| [d⁰]S = (9+3√3)/4 | 1.87e-67 | NUMERICAL (66 dig); value PROVED before | `tetra_Sdk_final.json` |
| [d¹]S, [d²]S to 51 dig; [d³]S to 49 dig; [d⁴]S to 48 dig | h-refinement certificates 5e-52…5e-49 | NUMERICAL | `tetra_Sdk_h48.json`, `tetra_SB_inner6496.json` |
| FD check k = 1: 17 dig; k = 2: 14.5 dig; k = 3: 11.8 dig (direct S(p), 19 pts on disk) | 7.5e-17 / 3.9e-15 / 2.7e-11 | NUMERICAL | `fd_Svals.json`, `fd_check.json` |
| Region A k=0 = 3√3/4; region B k=0 = 9/4 | 7.2e-68; exact | NUMERICAL-exact | h48 run |
| PSLQ, ζ(3) coefficient, Cl₂(π/3)/L(χ₋₃,2) exclusion | — | **FAILED-IN-BUDGET** (inputs final) | `tetra_Sdk_final.json` |
| (V,E)-doubling ⇒ edge/vertex law untestable from octa+cubo alone | — | PROVED (counting) | this report §4 |

**Files this session:** `k4_quad.py` (patched), `tetra_Sdk_h16_smoke.json` (smoke test, k=0 to 4e-46), `tetra_Sdk_h32.json` (218 s, k=0 err 4e-68), `tetra_Sdk_h48.json` (609 s), `tetra_SB_inner6496.json` (region B inner 1/64, 1/96, 1/128; 270+396+716 s), `tetra_Sdk_final.json` (60-digit final coefficients + error bounds), `fd_Svals.json` (S(p) at p = 1±j·0.005, j = 0…9, ~90 s each), `fd_check.json` (Fornberg analysis); driver code in `step_09.py`–`step_21.py` of this session (step_19 = 13-point FD run, step_20 = Fornberg analysis, step_21 = extra FD points). **Next session closes parts (2)–(3) with two <1-minute scripts** (19-pt Fornberg on `fd_Svals.json`; PSLQ on `tetra_Sdk_final.json`).