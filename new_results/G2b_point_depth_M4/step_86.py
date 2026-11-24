
G2B = "<project>/research/08_theory/G2b_point_depth_M4/"
open(G2B + "result.md", "w").write(r"""# G2b: M4 at the square, hexagon and c = 1/3 corners — the TRUE point-level depth test at weight 5

Work dir: `research/08_theory/G2b_point_depth_M4/`. Python: project `.venv`. Labels: **PROVED** (exact symbol algebra), **NUMERICAL**.

## Headline results

1. **The fixed point-level pipeline works; the control passes (PROVED).** G2's bug (pairing unexpanded exponent-vector words against basis-letter Lie words) is fixed by expanding with `cor6_sym.expand_sym` lifted to NF tensors (`step_01.py`). At c₀ = 1/√2 the true point-level weight-5 symbol of R₄'s weight-5 part (278 tensor words over the letters [c] = (−1/2,0), [1−c] = (−1/2,−1), [1+c] = (−1/2,1)) pairs to **(84, −38√2)** against `[X,[X,Y]]⊗[X,Y]` and `[Y,[X,Y]]⊗[X,Y]` with X = log 2, Y = log(1+√2) — exactly M4c's obstruction (normalisation factor +1; the pairing code itself reproduces (84, −38√2) from M4c's `symbol.json`, and the two symbol constructions were also checked to give the same tensor −[e₀]⁵ on the toy period J((0,0,0,0,2);1)). Shuffle products and Li₅-at-S-units pair to 0 (controls, 200 each).

2. **The hexagon corner is NOT classical at the symbol level (PROVED).** At c₀ = 1/2 (α = 2π/3; S-units {2,3}, X = log 2, Y = log 3):
   **obs(Ω₅(1/2)) = (41, −33/2) ≠ 0.**
   At c₀ = 1/3: **obs = (1838/3, −284) ≠ 0.** So M₄ needs depth two at the hexagon corner too — the letter specialisation does NOT kill the function-level obstruction of G2's step_12. (The naive specialised-function pairing at 1/2 is (11/2, −6): different from the true (41, −33/2), paralleling COR5's naive-vs-true discrepancy at weight 4.)

3. **Non-classical corners (sanity contrast, PROVED).** At c₀ = √3/2 (triangle; dimV = 3, basis {2, 3, −3(1+√3)}) and c₀ = (√5−1)/4 (decagon; dimV = 3, basis {2, 5, −3−√5}), all 24 pairings in Lie₃(V)⊗Λ²(V) are nonzero (`pair_nonclassical.json`).

4. **Depth-two correction at the hexagon, fixed by the cobracket (PROVED at the symbol level).** With the Goncharov cobracket formula (M4c's `delta_word`, ported; **exact check**: delta_word ≡ −δ⁻(5-symbol) mod shuffle products for every one of the 146 signed words, hence in T³/D³ for the total),
   Δ(R₄⁽⁵⁾(1/2)) = D ⊗ Li₂(1/3),  D ≡ (33/2)·Li₃(1/3) + 41·Li₃(−1/2)  (mod products; ζ(3)-part: see below),
   the Li₃-part being solved exactly in the quotient T³/D³. The minimal correction with CONVERGENT atoms (paper's convention Li_{3,2}(x,y) = Σ_{n>m>0} xⁿyᵐ/(n³m²); |x| < 1, |xy| < 1):
   **N₅ = (5/4)·Li₃,₂(1/4, 4/3) + (1/3)·Li₃,₂(8/9, 1)**,
   with Δ(N₅) = Δ(R₄⁽⁵⁾) in the quotient (coords (−41, −33/2)) and obs(N₅) = (41, −33/2) — both exact, re-run in `verify_g2b.py`. The ζ(3)-coefficient of D is invisible to the symbol; its consistency is witnessed numerically by the classical closed form of R₄⁽⁵⁾ − N₅ (search below).

5. **CONVENTION note (planner item 1).** Atom values are the n-outer iterated series Σₙ(xⁿ/n³)Σ_{m<n}yᵐ/m² = the paper's double series wherever absolutely convergent (|x| ≤ 1, |xy| ≤ 1); only such atoms were ever evaluated (166 values, 1300 digits). The first cobracket solve picked Li₃,₂(2, 1/3), which diverges in the paper's convention; **no value for it was ever computed or used** (step_04 only ever evaluated |x| ≤ 1 atoms). N₅ was re-solved over convergent atoms (`step_07.py`). Identity check: J((0,0,1/x,0,1/(xy));0→1) = **+**Li₃,₂(x,y), verified to 60 digits on six atoms (series both sides).

6. **Values (NUMERICAL).** M₄(2π/3) to 1300 digits from B4 (planner_check_b4's method; the He-sign bug in the first version — H_w(−c) needs H_w at −c, i.e. (−1)^nz(w)H_{−w}(c) — was found by the quadrature cross-check and fixed):
   M₄(2π/3) = 0.0021017551526063123975140574239779012704445906523150088079…
   **Cross-check vs G2's `values4_2pi_3.txt` quadrature: difference 4.87·10⁻⁹³** (= the quadrature's own precision). R₄⁽⁵⁾(1/2) = 0.0798379558673600680302842132879690163317943345890839233736… (1300 digits).

## Evidence table

| Check | Label | Result |
|---|---|---|
| control at 1/√2: pairing of Ω₅ == (84, −38√2) | PROVED (exact) | PASS |
| pairing code on M4c `symbol.json` | PROVED | (84, −38√2) |
| shuffle (200) and Li₅ (200 + 6 S-units) controls | PROVED | all 0 |
| hexagon c₀ = 1/2: obs(Ω₅) | PROVED | **(41, −33/2) ≠ 0** |
| c₀ = 1/3: obs(Ω₅) | PROVED | **(1838/3, −284) ≠ 0** |
| triangle, decagon: Lie₃⊗Λ² pairings | PROVED | 24/24 nonzero each |
| delta_word port == −δ⁻(5-sym) mod shuffles, all 146 words | PROVED | PASS |
| Δ(R₄⁽⁵⁾(1/2)) Li₃-part: (33/2)Li₃(1/3) + 41Li₃(−1/2) | PROVED | exact quotient solve |
| N₅ cobracket + pairing match (verify_g2b.py) | PROVED | PASS |
| M₄(2π/3) vs G2 quadrature | NUMERICAL | 4.87e-93 |
| J-word = +Li₃,₂ identity, 6 atoms | NUMERICAL | ≤ 1.8e-17 (N = 2500 series truncation) |
| classical closed form of R₄⁽⁵⁾(1/2) − N₅ | NUMERICAL | below (search status) |

## Files

Pipeline: `step_01.py` (fix + control), `step_02.py` (hexagon, c = 1/3), `step_03.py` (non-classical corners), `step_04.py` (1300-digit values; He-sign fixed), `step_05.py` (exact cobracket + first N₅ solve), `step_07.py` (convergent-atom N₅ re-solve), `verify_g2b.py` (all exact checks, ~1 s), `verify_g2.py` (B4 at two angles, ~3 s). Data: `pair_sqrt2.json`, `pair_half_third.json`, `pair_nonclassical.json`, `cobracket_hex.json`, `N5_convergent.json`, `values_hex.json` (208 values, 1300 digits), `progress.log`.
""")
print("result.md written")
