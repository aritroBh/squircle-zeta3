# THX1 status (auto-updated)

## DONE
- **MAIN THEOREM: Q5..Q9 all PROVED** with exact certificates (D~_{n,m} = zeta products, n=5..9 all m):
  - cert_D5_{0,1}, cert_D6_{0,1,2}, cert_D7_{0,1,2}, cert_D8_{0..3}, cert_D9_{0..3}.pkl saved.
  - Q_n(u) assembly matches Qpolys.pkl exactly (symbolic check).
  - verify_thx1.py: **24/24 PASS in 48 s** (exact Fraction verification of certificates for n<=8 + D9 m=0,
    mod-prime for D9 m=1..3; full exact checks of all D9 printed in main steps).
- Relation engine: DS+DUAL only (proved relations, no regularization). **Codim = DG bound at w<=9**
  (w6:13, w7:21, w8:34, w9:55) via flint rank (structs saved struct_w{5,6,7,8,9}.pkl).
  RUL1's RDS relations are FALSE (checked numerically, all 162 at w6) - excluded. DIST numerically OK but unused.
- Master GF derived + verified: D~(tau;m) = 2^{5/2} int (gd x)^{2m} H(x) [(1+e^{-x/tau})^tau - 1] dx.
- Individual I_{m,j} = j! m! G(0^{m+1}(-1)^j;1) reductions:
  - w<=5: already proved by RUL1 (their J-family). w6, w7: ALL certified (cert_I_{m}_{j}.pkl), PSLQ to 150 digits.
  - bases: w6: classical6 + zbar(5,-1); w7: classical7 + zbar(5,1,1) + zbar(4,2,1).
- w8 basis (34): classical8 (30) + zbar(6,1,1), zbar(5,2,1), zbar(5,1,1,1), zbar(4,2,1,1).
- w9 basis (55): classical9 (49) + zbar(7,1,1), zbar(6,2,1), zbar(5,3,1), zbar(5,1,1,1,1), zbar(4,2,1,1,1), X9f
  where X9f = G(0,1^6,0,-1;1) = -zeta(2,1,1,1,1,1,bar2). 247 {0,-1}-word representatives of its class found
  (w9_x9f_reps.pkl) for E1-computable numerics.
- Dixon p-adic lifting (thx1_dixon.py) is the workhorse for certificates; w9_Ainv.npy saved (Sf^{-1} reuse).

## DONE (cont.)
- Individual I_{m,j} reductions ALL CERTIFIED for w<=9:
  w<=5 (RUL1), w6 (5), w7 (6), w8 (7), w9 (8); certificates cert_I_{m}_{j}.pkl.
  Bases: w6 +zbar(5,-1); w7 +zbar(5,1,1)+zbar(4,2,1); w8 +zbar(6,1,1)+zbar(5,2,1)+zbar(5,1,1,1)+zbar(4,2,1,1);
  w9: classical9 + zbar(7,1,1)+zbar(6,2,1)+zbar(5,3,1)+zbar(5,1,1,1,1)+zbar(4,2,1,1,1)+X9fp,
  X9fp = G(-1,-1,-1,0,0,0,-1,0,-1;1) (E1-computable rep of the 55th direction, original X9f = -zeta(2,1,1,1,1,1,bar2)).
  NOTE: PSLQ candidates for I_{0,7},I_{1,6} were WRONG; corrected via modular route (free_reduce + basis solve).
  Numeric checks of every reduction to >= 107 digits.

## COMPLETE
- Assembly check: ALL OK (second proof route of Q5..Q9; non-zeta cancels exactly, zeta matches).
- V_n saturation computed: dim V_n = 3,3,4,4,5 for n=5..9; even n: V_n = span(e_m); odd n: + eta-singleton.
  Conjecture dim V_n = ceil(n/2) for all n  <=>  general Q_n reduction.
- Master GF proved: D~(tau;m) = 2^{5/2} int (gd x)^{2m} H(x) [(1+e^{-x/tau})^tau - 1] dx.
- verify_thx1.py: 43 PASS / 0 FAIL in 16 s.
- result.md + FINAL report written. 26 I-certificates + 16 D-certificates saved.
