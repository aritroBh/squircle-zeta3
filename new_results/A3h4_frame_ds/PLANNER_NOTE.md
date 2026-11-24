# Planner note for A3h4 (added while you run; please read)
- Do NOT use `A3f6_level16/evalseg.py`: it has a segment-order bug (it gives 0.66 for H(−1,−1,0), which should be ζ(3)).
  For level-16 numerics use `A3k_level16_funeq/mseval.py` (validated to 1e-121).
- `A3f6_level16/vals8_150digit.pickle` (level-8, weight 3) is fine. The planner checked H(−1,−1,0) = ζ(3) to 20 digits
  and the Kummer rows to 1e-139.
- There is no overlap with other jobs at weight 4 or for general-frame double shuffle. At weight 3, the pure inversion
  z ↦ 1/z family is covered by A3j (M1 session); skip it or just cross-check it.
