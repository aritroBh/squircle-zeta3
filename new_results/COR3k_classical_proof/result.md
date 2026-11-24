# COR3k: exact proof of the ten weight-3 HPL reductions and of B2′ — complete

**Deliverables (work dir `new_results/COR3k_classical_proof/`):** `verify_cor3.py` (+ its library `cor3_common.py`). End-to-end run: **34/34 items PROVED, wall time 74 s** (< 10 min), exit code 0, no stderr. Produced by steps `step_01.py`–`step_63.py` (exploration in 01–21; `cor3_common.py` built and patched in 22–51; component proofs tested in 52–58; `verify_cor3.py` written in 59/60; final full run in 62).

## Result: every identity in `classical_dict.json` is CORRECT. B2′ as typed in `result.md` is CORRECT.

### Method (exactly the prescribed algorithm; sympy never simplifies a polylog identity)
- Derivatives by chain rule; `polylog(1,a)` is rewritten **by hand** into `{log x, log(1−x), log(1+x), log 2}` (all factors > 0 on (0,1), hand factorizations like `log(2x/(1+x)) = log2 + log x − log(1+x)`).
- Vanishing = `sp.cancel() == 0` on every rational coefficient in `sp.Poly(..., domain=QQ(x))` over the atom symbols (π² tracked as an extra symbol).
- Constants at `x→0+`: hand-coded series with special values **substituted before comparison** (Li2: `0, π²/6, −π²/12, π²/12−log2²/2`; Li3: `0, ζ(3), −3ζ(3)/4, 7ζ(3)/8 − π²log2/12 + log2³/6`); the constant must be the **zero polynomial** in `Q[ζ(3), π, log 2]` — sufficient without any period-independence assumption.
- One algorithm repair was needed (documented): `d/dx E_w` is **not** weight‑1 (`R(x)·Li2(a)` differentiates to `R'·Li2(a) + logs`). Instead `E_w` is reduced by the **proved** weight‑2 atom identities to the basis `{Li2(x), Li2(−x), Li2((1+x)/2)}` + log monomials; every rational coefficient then cancels to 0 — a sufficient condition, no independence needed. (A hidden five-term relation `Li2((1−x)/(1+x)) = −Li2(x)+Li2(−x)+Li2((1+x)/2) −lx·l1m +lx·l1p +l1m·l1p −l1m·λ −l1p²/2 +λ²/2 +π²/12` was found as a nonzero coefficient, verified to 120 digits, then proved like the others.)

### Task 1 evidence (all EXACT)
| Item | Count | Result |
|---|---|---|
| Weight‑2 HPL identities `C_{a,b}=H_{a,b}`, b≠0 | 6 | derivative log‑coeffs all `0`; `C(0+)=0` zero‑polynomial → PROVED |
| Li2‑atom reductions (Euler/Landen/five‑term) + duplication `Li2(x)+Li2(−x)=½Li2(x²)` | 5+1 | derivative EXACT, value at 0+ EXACT → PROVED |
| Weight‑3 identities `H_w = C_w` | **10/10** | `d/dx C_w − f_{a1}C_{w'}` reduces to literally 0 terms; `lim_{x→0+} C_w = 0` → PROVED |
| Argument ranges | 8/8 | each arg in (−1,1) on (0,1), strictly monotone (grid min/max printed); principal branches real, chain rule valid; endpoint values {0, 1/2, ±1} only approached as limits |

### Task 2 evidence
- **Sign rule** `H_w(−x) = (−1)^{#nonzero} H_{−w}(x)`: one‑line proof (substitute `t→−t` in `H_{a,w}(x)=∫₀^x f_a H_w`; `f_0(−t)d(−t)=f_0(t)dt`, `f_a(−t)d(−t)=−f_{−a}(t)dt` for a=±1; induct from `H_{±1}(−x)=−H_{∓1}(x)`), numeric check 12/12 to 10⁻⁴⁰.
- `hpl_fit_k2_(0,0,0).json` ≡ typed B2 (`S(c)+S(−c)` via sign rule): all 12 coefficients in `Q[c,1/c]`, diff 0.
- **B2 → B2′**: substitution of the PROVED dictionary vs typed B2′ (minus `ζ(2)Φ(c)`, `Φ=½−q·artanh c`): **all 31 coefficient keys in Q(c) have diff = 0**, including `Z(c) = 7/4 + 2c`, `(c/3)λ³`, `−(c/6)π²λ`.

### NUMERICAL checks (not part of the proof)
| Check | Points | Residual | Digits |
|---|---|---|---|
| `H_w` quad of defining integral vs `C_w`, 10 words | x = 0.3, 0.7, 0.95 | 0 … 2.1·10⁻⁶⁵ | 65 dps (quad cross‑check vs fully nested: 1.8·10⁻⁶⁵) |
| typed B2′ vs stored `values_*.txt` M2 | π/2, 2π/3, π/3 | 1.1·10⁻¹²⁰, 5.2·10⁻¹²¹, 6.3·10⁻¹²¹ | 120 dps (100+ digits) |

### Corrections to the first attempt (`../COR3_classical_proof/`)
- Its η(3) constants were indeed exactly 0 (`η(3)=3ζ(3)/4`), but its weight‑2 candidates were genuinely wrong: its `C_{1,−1}` had `lim_{x→0+} = −log2²/2 ≠ 0`. Corrected (proved) forms used here:
  - `H_{1,−1}(x) = −Li2((1+x)/2) − log(1−x)log(1+x) + log2·log(1+x) − log2²/2 + π²/12`
  - `H_{−1,1}(x) =  Li2((1+x)/2) − log2·log(1+x) + log2²/2 − π²/12`

**Domain caveats:** all identities hold for `0 < x < 1`; no caveats on subintervals needed. **Labels:** everything above is PROVED except the two rows marked NUMERICAL.