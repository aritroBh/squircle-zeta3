import os, shutil
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd)
shutil.copy("<project>/research/08_theory/K1d_engine_fixed/esum_engine.py", "esum_engine.py")
md = """# DEC2 — the regular decagon's P'''(1): exact integrand, exact assembly (PROVED)

Status: Task 1 DONE (PROVED). Tasks 2-3 in progress.

## Task 1 — exact integrand (PROVED)
- Geometry derived from scratch over Q(sqrt5, rho), rho = sqrt(10-2 sqrt5) (step_01/02): normals
  n_i = (cos(i pi/5), sin(i pi/5)), edge V0 = (-h, y0), V1 = (-h, -y0) with h = (sqrt5-1)/4 = 1/(2 phi),
  y0 = (sqrt5-2)/rho; kink conditions <n2,V0> = <n3,V1> = 0 EXACT; sum_i a_i(t) = 1 EXACT (edge in dK1);
  term values (DEC1 order): a2 = t/(2 phi^2), a3 = (1-t)/(2 phi^2), a1 = (1+t/phi)/(2 phi^2),
  a4 = (1-t/phi^2)/(2 phi), a5 = 1/(2 phi) — all with coefficients in Q(sqrt5) (normal-form check);
  L^2 = (25-11 sqrt5)/10 EXACT (= paper's L5^2), h^2/L^2 = (5+2 sqrt5)/4, |y'|=L, y.y'/L^2=-(1-2t)/2.
- Length element (PROVED geometry): radial parametrization x(t) = lambda(t) y(t),
  lambda = F^{-1/(1+d)}, F = sum a_i^{1+d} (homogeneity), so
  |x'|^2 = F^{-2/(1+d)} ( |y'|^2 - 2 (y.y') Phi + |y|^2 Phi^2 ), Phi = S/F, S = sum a_i' a_i^d.
- Exact series engine (step_03): E3n = [d^3] sqrt(|x'|^2/L^2) computed in exact Q(sqrt5) arithmetic:
  55 log-monomials in the 5 term logs, t-degree <= 3.
- Monomial-by-monomial vs DEC1's reconstructed E3n_exact_pq.pkl: 55/55 keys, 440 rational
  coefficients, **0 mismatches** (step_04a).
- Independent check: direct numerical d-expansion of the radial-map speed at t = 0.3 (80 dps):
  |diff| = 1.03e-80 (step_04b).
- The only radical outside Q(sqrt5) is the overall edge length L = sqrt((25-11 sqrt5)/10)
  (the analog of the octagon's u), as expected.
"""
with open("result.md","w") as f:
    f.write(md)
with open("progress.log","a") as f:
    f.write("step_04c: result.md updated (Task 1 PROVED); esum_engine.py copied for exact J-tables\n")
print("ok")
