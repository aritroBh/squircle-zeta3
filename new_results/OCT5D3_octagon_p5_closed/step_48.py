import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
print("build log:", open("build_big.log").read())
print("cache:", [f for f in os.listdir(".") if f.startswith("gtab2_")])
# append note to result.md
md = open("result.md").read()
md += """
## Progress log
- Task 1 (exact symbol algebra): M4f basis built (70 atoms: 14 Li5 + 54 products + z5 + Lchi5;
  52 visible to the symbol). C = (-3600+2160sqrt2)Li32(b,1) + (11880-8280sqrt2)Li32(1/sqrt2,1):
  its (3,2) pair = (-1485+1035sqrt2, 900-540sqrt2) = the obstruction (PROVED, exact).
  S(W5)-S(C) lies EXACTLY in the span of the 52 visible M4f atoms (rank 30 over Q(sqrt2),
  kernel 22, residual tensor 0, obstruction 0) -> classical representative found (symfit_m4f.pkl).
- CONVENTION FIX vs OCT5D2: the value of Li32(x,1) consistent with the symbol convention
  gsym(0,0,1/x,0,1/(xy);1) is the GENUINE double polylog sum_{n>m} x^n y^m/(n^3 m^2).
  OCT5D2's atom_vals430 L32_b1/L32_h1 were off by a factor 1/x (off-by-one series); my values
  are smaller by exactly x (checked: ratio = b and 1/sqrt2 resp.). This bug likely contributed
  to OCT5D2's failed invisible fit.
- Local law m=4 (prop:order5, formulas cross-checked vs E5_edge_law/result.md):
  C_4^(5) = (33750-23805 sqrt2) cos(pi/8) = 78.2028433737847736987918987074...
  C_4^(2,3) = (-17370+12315 sqrt2) cos(pi/8) = 42.5354327315259866080186809336...
- gengine2p.py (parallel, 13 procs): dps=430 N=1440 rebuild 9.0s, residual 4.21e-432 (matches old).
  Big build dps=1650 N=5700 running.
"""
open("result.md","w").write(md)
print("result.md updated")
