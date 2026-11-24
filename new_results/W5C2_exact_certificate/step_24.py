import os, time
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
summary = """# W5C2 — exact certificate for W5 = N5 + R5_min3 (status: RUNNING)

Pipeline (mirrors A3h5 at weight 5):
- step_01: row trimming per sector (sparse echelon mod p1=16777213) -> S, C. [RUNNING]
- step_02: second-prime (p2=16777199) checks: rank(S)=|S|, det(R_{S,C})!=0, D in span(S).
- step_03: exact solve of y.R_{S,C}=D_C over Q per sector (Dixon/CRT), both RHS.
- step_04: exact Fraction verification on all sector columns; build certificate_{sec}_{comp}.json.
- verify_W5.py: standalone re-check from saved files (folds rebuilt exactly, per-sector equations,
  exact unfold identity on all 46656 columns).
"""
open('summary.md','w').write(summary)
if not os.path.exists('progress.log'):
    open('progress.log','w').write('')
print('summary.md written')