
with open('summary.md') as f:
    cur = f.read()
addition = """

## Kimi job (CHN2 Part B/C) — progress

### Part B cross-check: DONE (2026-10-29), full agreement
Exact symbolic expansion of (B+G0)/2 in u (log-Gamma series, Fraction arithmetic; polys in (x,y))
on Pi1 = (x,x,y,y) and Pi2 = (x,y,x,y), reduced to the CHN1 DS basis:
- W=2..11: EXACT polynomial identity with CHN1 forms.pkl restricted to each plane (all basis coordinates, both planes).
- W=12: EXACT agreement with engine.reduced_coeffs at 20 rational points per plane, all 12 basis coordinates
  (homogeneous degree-12 polys in (x,y) have 13 coeffs; 20 ratios decide identity). Irreducible coords
  z2^2 z(3,5), z2 z(3,7), z(2,1,1,8) vanish identically on both planes; z(3,9), z(5,7) appear in exact ratio
  z39 = -(14/9) z57 (the CHN1 weight-12 rescue), values match the engine.
- One self-correction during the work: my first Pi2 substitution had G0 wrong (u1+u3 = 2xu on Pi2, not (x+y)u);
  caught because Pi1 matched and Pi2 did not; after the fix everything matches. Planner theorem stands.
- Files: partB_planecoords.pkl; engine copies mzv_exact.py/engine.py/symengine.py; steps 13-18.

### Part C (P5): in progress
"""
with open('summary.md','w') as f:
    f.write(cur + addition)
print("summary.md updated, length", len(cur+addition))
