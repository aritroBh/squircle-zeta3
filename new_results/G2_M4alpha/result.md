G2 report — M₄ at every corner angle in HPLs (B4 found and **PROVED for all n**); depth-two obstruction at the function level (nonzero); point-level weight-5 computation unfinished (bug identified, see §3).

# G2: M₄ at every corner angle, in harmonic polylogarithms, and where it needs depth two

Work dir: `new_results/G2_M4alpha/`. Python: project `.venv`. Every number below was printed by a script in this folder. Labels: **PROVED** (exact algebra; all‑n identities by the Z‑sum method of `W1_corner_family/B_taylor_pi4/b3_allN_proof.py`), **EXACT through c^N**, **NUMERICAL** (digit counts given).

## 1. Task 1 — exact coefficients μ₄,ₙ (done)

**Theorem (coefficient content; PROVED, `step_05.py`).** μ₄,ₙ ∈ ℚ + ℚζ(2) + ℚζ(3) + ℚζ(4) + ℚζ(2)², nothing else: μ_{k,n} is a ℚ(n)‑combination of Bell polynomials of degree ≤ k in the cumulants κ_m ∼ ψ⁽ᵐ⁻¹⁾ (m ≤ k); for k = 4 the possible zeta monomials are exactly {1, ζ(2), ζ(3), ζ(4), ζ(2)²}. **ζ(5) and ζ(2)ζ(3) do NOT occur** (they would need degree ≥ 5).

- Reference table: G1's `flat_exact.py` (unmodified) at k = 4, N = 640 → `flat_exact_N640_K4.{json,out}`.
- Independent symbolic formula (`step_05.py`): μ₄,ₙ as an explicit 51‑monomial ℚ(n)‑polynomial in H⁽ʳ⁾_{n−1}, H⁽ʳ⁾_{2n−1}, r ≤ 4 (`mu4_formula.pkl`). **It agrees with `flat_exact.py` exactly for all n ≤ 640.**
- **Patterns (PROVED as rational‑function identities, `step_07.py`):**
  [ζ(4)]μ₄,ₙ = (n+3)/((2n+3)(4n²−1)); [ζ(2)²]μ₄,ₙ = −n/((2n+3)(4n²−1)); [ζ(3)]μ₄,ₙ = −[ζ(2)]μ₃,ₙ (the latter also as a Z‑sum identity in `step_11.py`). The μ‑formula's denominators are only {2n−1, 2n+1, 2n+3, n, n+1} — this forces the prefactor range e ∈ {−3,…,1} in B4.
- **Numerics (`step_06.py`, NUMERICAL):** Σ_{n≤400} μ₄,ₙ c^{2n} vs M4_800 (800‑digit squircle value): diff **2.18·10⁻¹⁰²** (truncation‑limited); vs `corner_series` at 2π/3: 1.53·10⁻⁹⁸, at 1.3 rad: 1.54·10⁻⁸⁵ (both quadrature‑limited). It reproduces the stated M₄ = 0.0060279931126091225045161082508776… to all its digits.

## 2. Task 2 — the HPL form B4 (done; PROVED for all n)

**Theorem B4.** For |c| < 1,
M₄(α) = R₄(c) + ζ(2)·Ψ₂(c) − ζ(3)·Ψ(c) + ζ(4)·F₁(c) + ζ(2)²·F₂(c), c = cos(α/2), with
- R₄ = `hpl_fit_k4_(0,0,0).json`: 109 terms, HPL weights 2–5 (3/8/28/70 terms at weights 2/3/4/5), prefactors c^e, e ∈ {−3,…,1}; unique solution (rank 606/606) of the 641‑equation even‑sector system;
- Ψ₂ = `hpl_fit_k4_(1,0,0).json`: 22 terms, weights 1–3 (unique, rank 79/79);
- Ψ = B3's ζ(2)‑part (G1 `hpl_fit_k3_(1,0,0).json`, 15 terms): the ζ(3)‑part of B4 is exactly **−ζ(3)·Ψ(c)**;
- F₁ = (7c/16 − 5/(8c) + 3/(16c³))·artanh c + 9/16 − 3/(16c²), F₂ = (−c/16 − 1/(8c) + 3/(16c³))·artanh c + 1/16 − 3/(16c²) (`b4_pi4_part.json`; artanh c = (H₁+H₋₁)/2).

**Status:**
- **PROVED for all n** (`step_11.py`): for each of the five parts separately, the even identity [c^{2n}](B4 − Σμ₄,ₙc^{2n}) ≡ 0 and the odd identity [c^{2n+1}]B4 ≡ 0 hold as exact identities of ℚ(n)‑combinations of alternating Z‑sums (canonical forms match identically in n: 29/5/2/1/1 Z‑words per part), and the monomial coefficients c^{−3}…c³ were verified directly. Both sides converge on |c| < 1, so B4 is an identity of analytic functions there. This upgrades B4 past G1's B3 (which was EXACT‑640 + numerics when found, proved later).
- **EXACT through c^1280** (`step_10.py`): all five combinations reproduce the exact target coefficients through c^1280; all odd coefficients vanish exactly.
- **NUMERICAL** (`step_10.py`, B4 evaluated from its own HPL Taylor series to order 7000 — probing coefficients far beyond the fitted c^1280): vs ~90‑digit `corner_series` M4 at six angles: 3.99·10⁻⁹² (π/6), 2.55·10⁻⁹² (π/4), 2.87·10⁻⁹² (π/3), 1.75·10⁻⁹³ (π/2), 4.87·10⁻⁹³ (2π/3), 1.56·10⁻⁹³ (1.3 rad); at c = 1/√2 vs M4_800: **9.64·10⁻¹²¹**.

## 3. Task 3 — depth structure (partially done; point level UNFINISHED)

**Function level (`step_12.py`, PROVED as an exact function computation):** the (3,2) part of the weight‑5 part of B4, modulo constants, **is NONZERO as a function of c**. Symbol in T⁵(V), V = {[c],[1−c],[1+c]} (constant letters such as [2] projected out); pairing in Lie₃(V)⊗Λ²V (dim 8×3) with w5_symbol.py's antisymmetrised convention; 18 of 24 pairings are nonzero Laurent polynomials, e.g.
[Y,[X,Y]]⊗X^Y = −3c⁻³+2c⁻¹+c; [Z,[X,Z]]⊗X^Z = +3c⁻³−2c⁻¹−c; [Y,[X,Z]]⊗X^Y = −(3/2)c⁻³+(5/2)c⁻¹−c; [X,[Y,Z]]⊗Y^Z = 2; [Y,[Y,Z]]⊗Y^Z = −4+4c (full matrix: `func_delta32.json`).
Controls: exhaustive — **all** shuffle products of total length 5 over 3 letters pair to zero under all 24 pairings (0 failures), plus 200 random products ×2 pairings in the script.
**Interpretation:** at the function level, M₄ already fails classicality — generically it needs depth two at *every* angle; this is not only a specialisation artifact of constant letters.

**Point level (UNFINISHED — bug identified, no results claimed).** The COR6 pipeline was adapted to weight 5 (`step_13.py`): certified S‑unit letter data printed and exact (at c₀ = 1/√2: field ℚ(√2), S = {2}, dimV = 2, basis {2, 1+√2} = M4c's alphabet; letter vectors [c] = (−1/2,0), [1−c] = (−1/2,−1), [1+c] = (−1/2,1); the weight‑5 symbol Ω₅(c₀) has 278 tensor words). **However the (3,2) pairing as coded returns empty: the symbol words must first be expanded into basis‑letter words (`cor6_sym.expand_sym`) before pairing with the Lie words — the current code pairs unexpanded exponent‑vector words.** Consequently the required control (reproduction of M4c's (84, −38√2)) and the answers at c₀ = 1/2 and c₀ = 1/3 were **never printed**; I make no claim about them. Note the tension to be resolved: the function‑level obstruction is nonzero, while the hexagon point could still be classical if the letter specialisation ([c]→−[2], [1−c]→−[2], [1+c]→[3]−[2]) kills it — this is exactly what the fixed `step_13.py` must decide.

## 4. Task 4 — hexagon closed form over MT(ℤ[1/6])

Not attempted: it was conditional on the point‑level classicality at c₀ = 1/2, which is undecided (see §3). No feasibility numbers were printed; nothing is claimed.

## 5. Evidence table

| Check | Label | Result |
|---|---|---|
| μ₄,ₙ content = {1,ζ2,ζ3,ζ4,ζ2²}; no ζ5, no ζ2ζ3 | PROVED (structural, `step_05.py`) | exact |
| symbolic formula == `flat_exact.py` (two independent codes) | PROVED (exact) | all n ≤ 640 agree |
| [ζ4], [ζ2²] formulas; [ζ3]μ₄ = −[ζ2]μ₃ | PROVED (rational‑function identities, `step_07.py`) | exact |
| Σ μ₄,ₙc^{2n} vs M4_800 / quadrature (2π/3, 1.3) | NUMERICAL | 2.18e‑102 / 1.53e‑98 / 1.54e‑85 |
| B4 all‑n identity, all five parts (even+odd Z‑sum identities + low orders) | **PROVED** (`step_11.py`) | all True |
| B4 == exact target through c^1280; odd coefficients vanish | EXACT‑1280 (`step_10.py`) | all True |
| B4 vs M4 at 6 angles (series order 7000 ≫ fitted c^1280) | NUMERICAL | ≤ 4.0e‑92 |
| B4(1/√2) vs M4_800 | NUMERICAL | 9.64e‑121 |
| function‑level (3,2) of B4's weight‑5 part | PROVED (exact; shuffle controls exhaustive) | **nonzero**, 18/24 pairings |
| point‑level (3,2) at 1/√2 vs (84,−38√2); at 1/2; at 1/3 | NOT DONE | bug: missing `expand_sym` before pairing; nothing printed |

## 6. Files

Pipeline (all in `new_results/G2_M4alpha/`): `step_05.py` (symbolic μ₄,ₙ + cross‑check vs `flat_exact.py`), `step_06.py` (series vs M4/quadrature; writes `values4_*.txt`), `step_07.py` (π⁴ part, PROVED), `g2_fit.py` (even‑sector fit machinery), `step_08.py` (ζ(3) check + ζ(2) fit), `step_09.py` (weight‑5 rational fit, 606 unknowns), `step_10.py` (assemble B4; EXACT‑1280 + numerics), `step_11.py` (all‑n proof), `step_12.py` (function‑level (3,2)), `step_13.py` (point level — **known bug**, see §3), `result.md`. Artifacts: `flat_exact_N640_K4{,_formula}.json`, `mu4_formula.pkl`, `mu4_zsum.pkl`, `hpl_fit_k4_(0,0,0).json`, `hpl_fit_k4_(1,0,0).json`, `b4_pi4_part.json`, `func_delta32.json`, `point_obs_{sqrt2,half,third}.json` (letter data only). Note: the environment auto‑saved my scratch cells as `step_01..step_56.py`; the numbered pipeline files above are the meaningful ones (they overwrote some scratch slots), and `step_11.py` imports `load_terms, pi4_terms` from `step_10.py`.

**verify_g2.py (< 5 min): NOT delivered.** The verification is currently distributed: run in order `step_05.py 640` (~1 min), `step_06.py` (~10 s), `step_07.py` (~20 s), `step_08.py` (~30 s), `step_09.py` (~100 s), `step_10.py` (~130 s), `step_11.py` (~3 s), `step_12.py` (~2 s) — about 5–6 min total; consolidating this into a single < 5‑minute `verify_g2.py` (e.g., N = 520 fits, skipping the reference‑table regeneration) and the one‑line `step_13.py` fix (`S = expand_sym(S, dim)` before pairing, likewise in `naive_obstruction`) are the two concrete remaining items, together with Task 4.