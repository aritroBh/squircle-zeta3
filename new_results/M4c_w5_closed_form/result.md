> **Status (planner, 2026-10-22 19:25 CDT).** Codex stopped at 19:14 (usage limit), after the search had finished and
> `w5_verify.py` had been written but not yet run. The planner ran it: PASS. Result:
> - W5 is **not classical**. Its exact symbol obstruction is (84, -38*sqrt2), which is nonzero.
> - W5 = R5 (`closed_form.md`, 70 terms) holds numerically to 1596 digits. Its depth-two part is the exact cobracket
>   correction N5 (4 Li_{3,2} terms).
> - The planner verification section is at the end of this file.

# M4c — W5: nonclassical point-level symbol (search in progress)

The requested classical expression has a nonzero necessary obstruction in the
point-level iterated-integral symbol. The replacement must include genuinely
nonclassical weight-five periods. An enlarged numerical basis is being tested.

Read the paper's current weight-four theorem and the M4a decomposition outputs.
The saved M4a value does not meet the requested precision. `w5_precision.py`
extends it using Arb ball quadrature, logarithmic endpoint substitutions, and
explicit endpoint-tail bounds, with roots rebuilt after setting precision.

`w5_precision.py` has completed both working precisions. See `precision.out` and
`W5_4200.json`, `W5_4400.json` for printed values and certified ball radii. The
saved strings agree beyond the requested decimal precision. The earlier failed
tail-bound run in `progress.log` is superseded by the finite-ball rerun.

`w5_symbol.py` computes the full constant-retaining symbol using only exact
arithmetic in Q(sqrt(2)). Write X=log(2), Y=log(1+sqrt(2)),
P=[X,Y], Q=[X,[X,Y]], and T=[Y,[X,Y]]. The brackets here mean commutators of
noncommutative tensor words. Pair the antisymmetrized deconcatenation with
Q tensor P and T tensor P. The output is `(84, -38 sqrt(2))` (unnormalized).
Every shuffle product and every classical Li5 has zero such pairing. The
script checks the product splits and the specified classical atoms separately.
The complete symbol and controls are in `symbol.json` and `symbol.out`.

This is an exact obstruction to a coaction-compatible classical identity. It
does not assert injectivity of the numerical period map, or rule out an
accidental numerical equality independently of that conjectural principle.
It does not rely on a moving-parameter function obstruction.

All writes remain in this directory. Computations run sequentially. Invocations
use `nice`, but this sandbox reports `setpriority: Operation not permitted`;
actual OS deprioritization cannot be claimed.

## Full point cobracket and an exact nonclassical correction

`w5_cobracket.py` computes the antisymmetric point coaction modulo products,
including the symbol-invisible primitive trilogarithms. The right weight-two
quotient is represented by Li2(b). Its left coefficient is

\[
\begin{aligned}
D={}&\left(-\frac{147}{2}+\frac{87\sqrt2}{16}\right)\zeta(3)
 -(374+76\sqrt2)\operatorname{Li}_3(b)\\
 &+(374-76\sqrt2)\operatorname{Li}_3(-b)
 +672\operatorname{Li}_3(1/\sqrt2).
\end{aligned}
\]

Define the depth-two convention
\(\operatorname{Li}_{3,2}(x,y)=\sum_{n>m>0}x^ny^m/(n^3m^2)\). Put

\[
\begin{aligned}
N_5={}&\left(779-\frac{695\sqrt2}{56}\right)\operatorname{Li}_{3,2}(b,1)
 -672\operatorname{Li}_{3,2}(1/\sqrt2,1)\\
 &+(-374+76\sqrt2)\operatorname{Li}_{3,2}(-b,1)\\
 &+\left(-405+\frac{4951\sqrt2}{56}\right)\operatorname{Li}_{3,2}(b,-1).
\end{aligned}
\]

The code proves exactly that the (3,2) Lie cobracket of \(W_5-N_5\) is zero.
These coefficients are exact coaction data; they do **not** yet establish a
closed form for the remaining classical part. This is a necessary obstruction
test, not an unconditional assertion that its kernel is classical.

The computation reuses the supplied weight-three logarithmic primitives and
their exact dilogarithm/trilogarithm reductions. Shuffle reduces each weight-three
GPL modulo products to a two-log integral. Every such projection is checked
against both primitive Lie-word symbol pairings; independent Li5 and shuffle
product cobracket controls also pass. `cobracket.out` records each reduction
identity and multiplier; `cobracket.json` holds the resulting coefficients.


## Planner verification (Claude, 2026-10-22 19:25 CDT)

**Verdict.**
- **Obstruction.** W5 has an exact symbol-level obstruction to classicality. So no identity W5 = (classical Li5's +
  products) holds at the motivic level. A purely numerical classical identity is ruled out only conditionally, by the
  usual period conjecture.
- **Closed form (numerical).** W5 = R5 = (Q(sqrt2)-combination of Li5, products, zeta(5)) + N5. Here N5 is four
  Li_{3,2} terms whose coefficients come exactly from the (3,2) cobracket.
- **Evidence.** R5 agrees with W5 to 1596 digits, with the exact symbol and the exact (3,2) cobracket matching. This is
  not an exact proof: weight-5 relation certificates are out of scope.

**What the planner ran or checked.**
1. `w5_verify.py` (Codex's; it had not run when Codex stopped). Output in `planner_w5_verify.out`. All PASS:
   - the expanded relation is identical at lattice precisions 1200 and 1500;
   - the full point symbol of R5 equals that of W5 exactly;
   - the four Li_{3,2} atoms agree with independent quadrature to 1e-92;
   - R5 - W5(5400 bits) = -2.1e-1596;
   - the full (3,2) Lie cobracket of R5 - W5 is 0.
2. `planner/planner_symbol_pairing.py` (independent code; reads only `symbol.json`):
   - it recomputes the antisymmetrised (3,2) deconcatenation paired with [X,[X,Y]]⊗[X,Y] and [Y,[X,Y]]⊗[X,Y];
   - it gets **84** and **-38*sqrt2**, matching Codex;
   - it checks that 200 random shuffle products and 200 random Li5 symbols pair to 0.
3. `planner/planner_eval_R5.py` (independent mpmath code; no M4c imports):
   - Li_k uses mpmath.polylog, and Li_{s,t} uses a direct double series with the stated convention.
   - R5 - W5(M4a W5_600.txt; Kimi's separate quadrature) = 2.8e-601. That is the precision limit of W5_600.
   - R5 - W5(M4c 5400 bits) = 9e-837 at 820 digits of working precision.
4. `planner/planner_assemble_M4.py`: **full closed form of M4**, M4 = R5 + R''_nice + W3^(4) + T.
   - R''_nice is M4b2's proved form. W3^(4) = W3' - zeta(3)/8 - pi^2*lambda/24 - lambda^3/12, with W3' from thm:M3.
     T = (pi^2 + 6*lambda^2)/24.
   - Result: closed form - M4_800 (M4a, full-integrand quadrature) = **1.2e-803**, the precision limit of M4_800.

**Interpretation.**
- The first nonzero (3,2) cobracket component is the point where M_k stops being classical.
- M0 to M3 are classical: M3 needs Li4 at seven arguments.
- M4 is the first squircle coefficient that needs a genuine depth-two polylogarithm, Li_{3,2}.
- The depth-two coefficients in R5 match the exact cobracket correction N5 term for term. The LLL step was free to
  move them inside the symbol kernel, so this is an independent consistency check of the cobracket computation.

**Open.**
- R5 is ugly: integer-relation height about 1.4e15, denominators up to 2.3e7. This comes from the Li5 arguments the
  greedy step picked, e.g. -16+12*sqrt2 and 1/2-3*sqrt2/8. A nicer form, with Li5 at W4's seven arguments plus a few
  more and the same N5, is a mechanical follow-up.
- An exact weight-5 certificate (as for W4) is open.
