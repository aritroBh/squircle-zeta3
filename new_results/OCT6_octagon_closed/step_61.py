md = open('result.md').read()
old = """LLL at 745 digits then isolated each A-atom; the two reductions sum to the proved Landen
identity (cross-check 1.0e-743). Λ = Li₂₁(x1, −√2/2) and swap = Li₂₁(x2, −√2) reduce likewise
(LLL, 745 digits; Λ cross-validated by the stuffle Li₂(x1)Li₁(−√2/2) = Λ + Li₁₂(−√2/2, x1) + Li₃(x2)
with Li₁₂(−√2/2, x1) from an independent double series: residual 6.7e-744)."""
new = """LLL at 745 digits isolated each A-atom (R1, R2 below). Both are now PROVED analytically:
three function identities, each derived by differentiation (the derivative identities collapse
by dilog Landen/Euler/inversion) and pinned at explicit points, hold simultaneously —
  (L2, Landen sum)      A(x) + A(x/(x−1)) = z2·log(1−x) + 2log(1−x)Li₂(1−x) + ½log²(1−x)·logx − 3Li₃(1−x) + 3z3
  (DIFF)                A(x) − A(1−x)     = 2Li₃(x) − 2Li₃(1−x) + log(1−x)Li₂(1−x) − logx·Li₂(x) − ½log²x·log(1−x) + ½log²(1−x)·logx
  (L1, inversion)       A(−y) + A(−1/y)   = −z2·logy − (log³y)/3 + Li₃(−y) + 2Li₃(−1/y) − logy·Li₂(−y) + z3
verified at 3–5 generic points to 400–745 digits each. The orbit system at {x1, x2, √2/2, −(1+√2)}
is a solvable 4×4: A(x1) = [R_L2(x2) − DIFF(x2) − R_L2(√2/2) + R_L1(κ1)]/2, A(x2) = R_L2(x2) − A(x1);
these match the direct series to 7.3e-746 / 1.1e-746 and the LLL forms R1/R2 to 742+ digits.
Λ = Li₂₁(x1, −√2/2) and swap = Li₂₁(x2, −√2) reduce likewise (R3, R4: LLL-discovered, verified
745 digits; Λ cross-validated by the stuffle Li₂(x1)Li₁(−√2/2) = Λ + Li₁₂(−√2/2, x1) + Li₃(x2)
with Li₁₂(−√2/2, x1) from an independent double series: residual 6.7e-744)."""
assert old in md
md = md.replace(old, new)
# also update the R1/R2 resid line note
md = md.replace("(resid 1.3e-742)", "(PROVED; resid 1.3e-742)")
md = md.replace("(resid 1.4e-743)", "(PROVED; resid 1.4e-743)")
open('result.md','w').write(md)
print("result.md updated", len(md))
import os
steps = sorted(f for f in os.listdir('.') if f.startswith('step_') and f.endswith('.py'))
print(len(steps), "step files:", steps[:8], "...", steps[-4:])
