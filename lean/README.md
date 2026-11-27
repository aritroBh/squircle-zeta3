# Lean checks for "Investigating ζ(3)"

Lean 4 (v4.36.0) with Mathlib (v4.36.0). This folder machine-checks 565 statements taken from the paper (562 distinct: three in `Wave7/V04` repeat `Wave3/W04`), in 68
files. Among them, in full: `thm:M1` and `thm:M2` starting from the curve, `thm:rhombus`, `thm:zero`, `thm:fn`,
`thm:fnangle`, `thm:fourier`, `thm:xnyn`, `thm:poly`, `thm:P2`, `thm:one`, `thm:flat`, `thm:pip`, `thm:decay`, `thm:varpi`,
`thm:zladder`, `prop:degree`, and `thm:pw` (squircle, every rhombus, every dimension). It does not check every theorem;
`LEAN_CHECKLIST.md` says what is in and what is out.

- `LEAN_CHECKLIST.md`: every target, its status, and what Lean cannot do here.
- `Zeta3/Wave1/` to `Zeta3/Wave10/`: one file per target (`Wave2/D00_Dilog.lean` is the dilogarithm library). Each file quotes the paper, states the claims as `def stmt_k : Prop`
  between the `BEGIN STATEMENTS` and `END STATEMENTS` markers, and proves them below.
- `audit/`: the statement blocks as accepted by an independent audit against the paper. Proofs may not change them.
- `tools/gate.py`: the mechanical gate. `python3 tools/gate.py check Zeta3/Wave1/<file>` passes only if the
  statements are the audited ones, the file compiles, there is no `sorry`, and the only axioms are
  `propext`, `Classical.choice` and `Quot.sound`.

To rebuild: install `elan`, then in this folder run `lake exe cache get` and `lake build`.
