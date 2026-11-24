import os
base = '<project>/research/08_theory/'
summary = r"""# PG6: 14-gon (2,1) invariant — running summary

## Task A — P'''(1) discrepancy: RESOLVED
- CORRECT value: **P'''(1) = -5.5080836980534783** (the exact-integrand chain 84*L*int E3).
  Confirmed to ~24 digits by three independent computations:
  - high-precision int E3n = -0.645538287054672257864..., x84*L = -5.508083698053478284
  - Fourier filter on the FULL (sliver-free) chamber, rho=0.015: -5.508083698053478284 (24-digit agree)
  - Richardson-extrapolated finite differences: -5.508083697 (O(h^2) trunc)
- The polar-perimeter GOLD value -5.5080732424212877 was WRONG.
- BUG (located): step_06_gold14.py / step_06b_gold14.py integrate perim over
  [pi/2+1e-9, pi/2+pi/7-1e-9], removing 1e-9 slivers at the VERTEX endpoints. Near the vertex the
  p-ball integrand has a boundary layer: |cos(th-pi/2)|^p with th-pi/2 in [0,1e-9] multiplies each
  p-derivative by ~ln(1e-9) = -20.7, so the missing sliver contributes err(p) with err'''(1) ~ +1.05e-5.
  Being a fixed systematic error (same for every p-sample), the Fourier filter at rho=0.02 AND 0.01 and
  the FD stencil all extract the SAME err'''(1), hence appear self-consistent (7e-18) yet are all wrong.
  Verified: gold chamber(slivers) reproduces -5.5080732424; full/centered chambers give the truth.
- CONTROL (12-gon, known P'''(1)=-5.2032404295176662): corrected pipeline (chamber split at the kink
  theta=0, NO slivers) gives P(1)=12sqrt(26-15sqrt3) to 5e-70, Fourier P3 to 20 digits, FD to 2e-11.
  Note: 12-gon chamber center theta=0 has a vanishing |sin th| term (normal at pi/2) -> kink, needs split.
  The 14-gon has NO normal at pi/2, so its chamber center is smooth (no split needed).

## Task B — delta(N3): in progress
"""
open('summary.md','w').write(summary)
print("summary.md written")
# Explore DOD1 files (the engine)
for d in ['DOD1_dodecagon','DOD5_dodecagon_termwise','PGEN1R_referee_alphabet']:
    dd = os.path.join(base,d)
    print("="*30, d)
    for f in sorted(os.listdir(dd)):
        if f.endswith('.py') or f.endswith('.md'):
            print(f"  {os.path.getsize(os.path.join(dd,f)):>8}  {f}")
