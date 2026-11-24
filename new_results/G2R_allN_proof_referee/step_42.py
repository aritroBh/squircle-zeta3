CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
md = open(CWD + "/result.md").read()
old = """Corroboration reproduced in this folder: `step_10.py` exact check **EXACT through c^1280, all five
parts: True** (30 s); `zsum.py` self-tests True; planner's independent
`planner_check_b4.py`: |B4 − M4_800| = **8.6732e-131** at c = 1/√2 (dps 130), 4.8651e-93 at 2π/3,
2.8715e-92 at π/3; step_07's [ζ(4)]μ = (n+3)/((2n+3)(4n²−1)) and [ζ(2)²]μ = −n/((2n+3)(4n²−1))
re-derived from the pickle: both True."""
new = """Corroboration reproduced in this folder: `step_10.py` in full (147 s): **EXACT through c^1280, all
five parts: True**, and its series-order-7000 numerics match the claimed values: 3.99e-92 (π/6),
2.55e-92 (π/4), 2.87e-92 (π/3), 1.75e-93 (π/2), 4.87e-93 (2π/3), 1.56e-93 (1.3 rad),
|B4 − M4_800| = **9.64e-121** at c = 1/√2 (120 dps); `zsum.py` self-tests True; planner's
independent `planner_check_b4.py` (own mpmath evaluator, dps 130): |B4 − M4_800| = **8.6732e-131**
at c = 1/√2, 4.8651e-93 at 2π/3, 2.8715e-92 at π/3; step_07's
[ζ(4)]μ = (n+3)/((2n+3)(4n²−1)) and [ζ(2)²]μ = −n/((2n+3)(4n²−1))
re-derived from the pickle: both True."""
assert old in md
md = md.replace(old, new)
old2 = "- `run_base/` (pristine pipeline copy) with `run_base.out` (step_11) and `run_base_step10.out`\n  (EXACT-1280); `mut_a/`…`mut_d/` with `mut_a.out`…`mut_d.out` (the four mutants)."
new2 = "- `run_base/` (pristine pipeline copy) with `run_base.out` (step_11) and `run_base_step10.out`\n  (full step_10: EXACT-1280 + 6-angle numerics); `mut_a/`…`mut_d/` with `mut_a.out`…`mut_d.out`\n  (the four mutants)."
assert old2 in md
md = md.replace(old2, new2)
open(CWD + "/result.md", "w").write(md)
with open(CWD + "/progress.log", "a") as fh:
    fh.write("step_result_final: result.md updated with full step_10 numerics (9.64e-121 at the square); referee complete, verdict PROVED\n")
print("result.md updated; final listing:")
import os
for f in sorted(os.listdir(CWD)):
    p = os.path.join(CWD, f)
    if os.path.isfile(p):
        print(f"  {f}  {os.path.getsize(p)}")
    else:
        print(f"  {f}/  dir")
