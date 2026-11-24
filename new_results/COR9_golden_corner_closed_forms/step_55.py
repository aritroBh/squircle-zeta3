WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
r = '''# COR9 -- number-level closed forms of M3 at the three remaining classical corners + larger census

Work dir `research/08_theory/COR9_golden_corner_closed_forms/`. Labels: **PROVED** = exact algebra;
**NUMERICAL** = lattice relations identical at two precisions (dps 700 and 1000), residual at working
precision (dps 1120). Values: HPL = G1's proved B3 term lists (M3 = R3 + zeta(2) Psi - zeta(3) Phi)
by exact series at dps 1150; quad = independent tanh-sinh quadrature (COR7's precision-boosted
integrand). Machinery: `cor9_lib.py` (reuses cor7_lib/cor7_pt/cor6; adds CF-based fundamental units,
diop_DN norm-l elements, class numbers by reduced-form cycles -- validated on 62 fields).
Notation: b = sqrt2-1, L2 = log(1+sqrt2) (Task 1) / L = log phi (Task 2), lam = log 2.

## Task 1: c0 = (sqrt2-1)^2 = 3 - 2 sqrt2 (alpha = 2.79674065816409510609187082939 rad)

- **M3 = 0.001180612468597228774341058761863427649592...** (1145 digits, `b2sq_M3_hpl_1150.txt`).
  Two methods: |quad(dps 100) - HPL(1150)| = 4.33e-103 (>= 60 digits required). NUMERICAL.
- Weight split (COR5 convention) at dps 1150: W4 = 1.38878210588956677292055787703,
  W3 = -2.19524014103291099661881673522, W2 = 0.80763864761194145247259991695.
- **W2** = -Li2(c) - Li2(-c) + pi^2/12 = -(1/2) Li2(b^4) + pi^2/12 (duplication). Residual 0 @1120. PROVED (from the HPL formula + duplication).
- **W3** (COR7 12-element Q(sqrt2) basis x {1,sqrt2}; identical dps 700/1000; residual 2.28e-1116 @1120; a0 = 48). NUMERICAL:
  W3 = (-273/2 - 1655/16 s2) zeta3 + (-1166 - 898 s2) Li3(b) + (562 + 422 s2) Li3(-b) + (1544 + 1180 s2) Li3(1/s2)
     + (-314 - 236 s2) L2 Li2(b) + (392 + 295 s2) lam Li2(b) + (-747/8 - 287/4 s2) pi^2 L2 + (249/4 + 767/16 s2) pi^2 lam
     + (-193/6 - 295/12 s2) lam^3 + (-193 - 295/2 s2) lam^2 L2 + (272 + 413/2 s2) lam L2^2 + (56/3 + 47/3 s2) L2^3.
- **W4** (thm:W4R 29-element basis x {1,sqrt2}; 59-dim lattice; identical dps 700/1000; residual 1.65e-1115 @1120; a0 = 5760; height 2753280). NUMERICAL:
  Li4 part: (-478 - 153 s2) Li4(b) + (-470 - 153 s2) Li4(-b) + (83/2 + 59/4 s2) Li4(-b^2) + (4 - 4 s2) Li4(1/s2)
     + (103/2 + 41/4 s2) Li4(1/2) + (16 + 14 s2) Li4(1-1/s2) + (116 + 30 s2) Li4(2-s2);
  Li3/log part: (12 + 8 s2) L2 Li3(b) + (-8 - 4 s2) L2 Li3(1/s2) + (-9/4 - 3/2 s2) L2 zeta3
     + (-14 - 10 s2) lam Li3(b) + (8 + 4 s2) lam Li3(1/s2) + (17/4 + 21/8 s2) lam zeta3;
  weight-2 products: (-15 - 10 s2) Li2(b)^2 + (-6 - 4 s2) Li2(b) L2^2 + (-8 - 5 s2) Li2(b) lam L2 + (5 + 3 s2) Li2(b) lam^2 + (13/6 + 17/12 s2) pi^2 Li2(b);
  pure products: (-361/320 - 1681/5760 s2) pi^4 + (7/3 + 5/4 s2) pi^2 L2^2 + (25/6 + 11/24 s2) pi^2 lam L2 + (-43/12 - 3/4 s2) pi^2 lam^2
     + (-1/4 - 25/24 s2) L2^4 + (-52/3 - 9/2 s2) L2^3 lam + (61/4 + 35/8 s2) L2^2 lam^2 + (-9/2 - 5/8 s2) L2 lam^3 + (235/96 + 13/64 s2) lam^4.
- **Lie class (PROVED):** the Li4 part of the numerical W4 has Q^4-class components (84, -100, 128) and
  (-4, -16, 48) = exact class(Omega((sqrt2-1)^2)) (cor7_pt/cor6 symbol algebra) = COR7's printed class.
- Assembly: |M3 - (W4cf+W3cf+W2cf)| = 8.54e-1118 @1120.

## Task 2: golden corners, field Q(sqrt5), S = {2} (2 inert; S2-units = <-1, phi, 2>)

Period bases over Z[phi, 1/2], coefficients in Q(sqrt5) (1 and sqrt5 parts). S-unit arguments x
(x and 1-x both {2}-units) enumerated exactly: 27 (K = 9 box; COR7's 25 with K = 8). Li4 classes of the
27 args span Q^3 (rank 3). Exact relations found and removed before LLL:
- weight 3 (relation rank 3 among 8 candidate values): Li3(1/phi^2) = 4 Li3(1/phi) + 4 Li3(-1/phi) (duplication);
  6 zeta3 = 30(Li3(1/phi)+Li3(-1/phi)) + pi^2 L - 5 L^3; and a third relation (r_c - r_d in step_39).
  Independent weight-3 core: {zeta3, Li3(1/phi), Li3(1/phi^3), Li3(1/(2phi)), Li3(2/phi^2)} + products.
- weight 4 (relation rank 4 among the 12 candidate Li4): Li4(1/phi^6) = 8 Li4(1/phi^3) + 8 Li4(-1/phi^3)
  (duplication x8) plus 3 more (heights 288, 567, 642; step_44). Independent Li4 set (classes still span Q^3):
  {-sqrt5 phi, 1/(2 phi^2), 1/phi^2, 1/2, 2/phi^2, phi/2, 1/phi^6, -1/phi^3}.
- Li2 is elementary: Li2(1/phi) = pi^2/10 - L^2, Li2(1/phi^2) = pi^2/15 - L^2, Li2(-1/phi) = -pi^2/15 + L^2/2
  (residuals ~1e-302). Li2(1/phi^3) is NOT elementary (no small relation at dps 250) -> its products are in the basis.

### c0 = 1/phi (alpha = 1.809113788604762728254634 rad)

- **M3 = 0.01695853531118527475323767966725547618373...** (1145 digits). |quad(dps 100) - HPL| = 6.7e-103. NUMERICAL.
- W4 = 0.9523025492009508271656064, W3 = -1.544606644232830679602558, W2 = 0.6092626303430651271901891 (dps 1150).
- **W2** = -Li2(1/phi) - Li2(-1/phi) + pi^2/12 = **pi^2/20 + L^2/2** (elementary). Residual 0 @1120. PROVED.
- **W3** (13-element basis; identical dps 700/1000; residual 6.93e-1119 @1120; a0 = 240). NUMERICAL:
  W3 = (4/5 - s5/10) zeta3 - Li3(1/phi) + (-41/240 - s5/16) pi^2 L + (-7/24 - 5 s5/24) L^3.
- **W4** (32-element basis; 65-dim lattice; identical dps 700/1000; residual 2.18e-1115 @1120; a0 = 64800; height 1664550). NUMERICAL:
  Li4 part: (-5/2 + 11/6 s5) Li4(1/(2 phi^2)) + (113/16 - 187/96 s5) Li4(1/phi^2) + (3/10 + 19/30 s5) Li4(1/2)
     + (-7/2 + 3/2 s5) Li4(2/phi^2) + (-2 - 2/3 s5) Li4(phi/2) + (-29/72 + 1/12 s5) Li4(1/phi^6) + (35/18 + 7/6 s5) Li4(-1/phi^3);
  rest: (339/40 + 113/40 s5) L zeta3 + (-1 - 2 s5) L Li3(1/phi) + (-5 + 2/3 s5) L Li3(1/phi^3) + (-5 - s5) L Li3(2/phi^2)
     + (-7/30 + 1/60 s5) pi^2 Li2(1/phi^3) + Li2(1/phi^3)^2 + (25/4 + 3/4 s5) Li2(1/phi^3) L^2 - 6 Li2(1/phi^3) L lam
     + (64/2025 - 1/5400 s5) pi^4 + (-77/40 - 67/180 s5) pi^2 L^2 + (1/5 + 31/180 s5) pi^2 L lam + (3/10 - 1/30 s5) pi^2 lam^2
     + (411/16 + 2041/288 s5) L^4 + (-24 - 143/18 s5) L^3 lam + (6 + 5/3 s5) L^2 lam^2 + (-1/2 + 13/18 s5) L lam^3 + (-7/40 + 3/40 s5) lam^4.
- **Lie class (PROVED):** Li4 part = ((111/2, 0, -12), (-51/2, 0, 12)) = class(Omega(1/phi)) = COR7's printed class.
- Assembly: |M3 - (W4cf+W3cf+W2cf)| = 2.13e-1120 @1120.

### c0 = sqrt5 - 2 = 1/phi^3 (alpha = 2.664957729970061020416 rad)

- **M3 = 0.002253617416574831375768762960540799083009...** (1145 digits). |quad(dps 100) - HPL| = 3.6e-103. NUMERICAL.
- W4 = 1.35757721488607082836685, W3 = -2.149528455790214148331482, W2 = 0.7942048583207181513404015 (dps 1150).
- **W2** = -Li2(1/phi^3) - Li2(-1/phi^3) + pi^2/12 = -(1/2) Li2(1/phi^6) + pi^2/12. Residual 0 @1120. PROVED.
- **W3** (13-element basis; identical dps 700/1000; residual 7.55e-1118 @1120; a0 = 120). NUMERICAL:
  W3 = (-2029/40 - 1061/40 s5) zeta3 + (-19 - 17/2 s5) Li3(1/phi) + (35/3 + 5 s5) Li3(1/phi^3) + (62 + 32 s5) Li3(2/phi^2)
     + (151/10 + 116/15 s5) pi^2 L + (-47/6 - 4 s5) pi^2 lam + (-48 - 24 s5) Li2(1/phi^3) L + (32 + 16 s5) Li2(1/phi^3) lam
     + (-3905/12 - 2003/12 s5) L^3 + (699/2 + 180 s5) L^2 lam + (-93 - 48 s5) L lam^2.
- **W4** (32-element basis; identical dps 700/1000; residual 2.49e-1116 @1120; a0 = 32400; height 693000). NUMERICAL:
  Li4 part: (1/6 + 1/2 s5) Li4(1/(2 phi^2)) + (-163/48 + 13/16 s5) Li4(1/phi^2) + (29/30 - 11/10 s5) Li4(1/2)
     + (11/2 + 1/2 s5) Li4(2/phi^2) + (11/3 + s5) Li4(phi/2) + (-1/8 + 1/72 s5) Li4(1/phi^6) + (1 + 7/9 s5) Li4(-1/phi^3);
  rest: (52/5 + 18/5 s5) L zeta3 + (-27/5 - 8/5 s5) lam zeta3 + (-9 - 4 s5) L Li3(1/phi) + (4 + 2 s5) lam Li3(1/phi)
     + (-191/2700 - 47/8100 s5) pi^4 + (-155/72 - 13/24 s5) pi^2 L^2 + (113/45 + 1/3 s5) pi^2 L lam + (-2/3 - 1/30 s5) pi^2 lam^2
     + (1939/144 + 59/48 s5) L^4 + (-385/18 + 1/6 s5) L^3 lam + (133/12 - 3/4 s5) L^2 lam^2 + (-43/18 + 1/6 s5) L lam^3 + (17/60 + 1/60 s5) lam^4.
- **Lie class (PROVED):** Li4 part = ((84, -21, 6), (-18, 3, 0)) = class(Omega(1/phi^3)) = COR7's printed class.
- Assembly: |M3 - (W4cf+W3cf+W2cf)| = 3.26e-1119 @1120.

## Task 3: census D <= 200, heights |a|, |b|, e <= 40  (in progress)
'''
open(WD+'result.md','w').write(r)
print("result.md written,", len(r), "chars")
open(WD+'progress.log','a').write(f"result.md updated with Tasks 1-2.\n")
