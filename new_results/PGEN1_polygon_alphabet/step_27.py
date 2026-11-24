import time
open("progress.log","a").write(f"{time.strftime('%Y-%m-%d %H:%M:%S')} step_01..08: read paper subsection + OCT9c/OCT9d/DOD5/DEC1/OCT6 result.md; located exact data (short_form.pkl, TW3.pkl, tables in result.md's)\n")
open("result.md","w").write("""# PGEN1: the trilogarithm alphabet of a general polygon, and the cobracket of the far kink lines

(work in progress; updated as steps complete)

## Plan
1. Alphabet rule: check on hexagon (thm:hexagon), generic rational octagon (per edge, OCT9d data), regular octagon, decagon, 12-gon; symbol-level proof for s <= 4 far kink lines.
2. delta(N_3) = sum c_i {z_i}_2 (x) z_i computed via wedge (z^(1-z)) (x) z in ^2(F^x) (x) F^x, for: regular octagon, generic octagon (edge by edge and total), decagon, 12-gon; controls first.
3. Optional: single-valued L3 on OCT9c vs OCT9d bases.

## Data sources (read-only)
- OCT9d short_form.pkl: exact Li3 coefficients per edge (PROVED symbol identities M x = TW3[e]): G8 = [-9/16, 4/9, 1/5, 1/3, -5/3, -1/5, 5/4, 2/5] for edges 1-3, B0 (20 args) for edge 0.
- OCT6 result.md: regular octagon Li3 coefficients in S (P'''(1) = 48 u S).
- DEC1 result.md: decagon Li3(1/phi), Li3(-1/phi) coefficients.
- DOD5 result.md: 12-gon Li3 coefficients at 8 points.
- kappa values per octagon edge: K12 = {0:(3/2,-7/12), 1:(1/3,1/9), 2:(5/3,-5/9), 3:(1/4,1/24)} (oct9c_common.py).
""")
print("wrote result.md and progress.log")
