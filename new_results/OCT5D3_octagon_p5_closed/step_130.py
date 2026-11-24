import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
md = open("result.md").read()
md += """

## FINAL PROVED SYMBOL CHECK
S(closed form of Wt5) = S(W5oct) EXACTLY (difference tensor 0), and its (3,2) obstruction
= (-1485+1035sqrt2, 900-540sqrt2). So the classical+Li32 part is symbol-correct (PROVED);
the 21 invisible products have zero symbol by construction. Exactness of the full identity
is NUMERICAL (modulo the period conjecture), verified to ~1490 digits at two precisions.

## Evidence table (final; all numbers printed from code)
| # | Check | Result | Residual | Digits | Producer |
|---|---|---|---|---:|---|
| 1 | obstruction S(W5oct) | PROVED | exact | - | verify_oct5d3.py A |
| 2 | C's pair = obstruction | PROVED | exact | - | verify A |
| 3 | S(W5)-S(C) in span of 52 M4f atoms, residual 0 | PROVED | exact (rank 30, kernel 22) | - | verify A |
| 4 | Li32 symbol kernel dim 1 (K32) | PROVED | exact | - | verify A |
| 5 | S(closed form Wt5) = S(W5oct) | PROVED | exact, obstruction matches | - | wt5_primitive step |
| 6 | G-table dps 1650 N 5700 | NUMERICAL | 3.54e-1652 | 1652 | gengine2p, build_big.log |
| 7 | Wt0 = 0 | PROVED (QS2) | exact | - | verify B |
| 8 | sum(Wt_k) vs p5_800 P5/u | NUMERICAL | 8.24e-800 | 800 | assembly step |
| 9 | Wt1 exact form | PROVED + num | 0.0 (LLL dps 800) | 800 | Wt1 step, verify B |
| 10 | Wt2, Wt3, Wt4 closed forms, identical 800(900)/1500 | NUMERICAL | 1.4e-1497 / 1.9e-1497 / 1.5e-1495 | 1497 | q2_lll steps; verify B (1.3e-349) |
| 11 | Wt5 closure rel1500, direct eval 900/1100/1500 | NUMERICAL | 1.17e-894 / 1.26e-1090 / 1.4e-1494 | 1490 | closure steps; verify B (1.1e-351) |
| 12 | rel1500 - rel1100 in verified atom-kernel | exact rank check | rank 26 = 26 | - | kernel-membership step |
| 13 | closed form reproduces P^(5)(1), P5/u | NUMERICAL | 4.88e-794 / 2.64e-794 | 794 | assembly step; verify C |
| 14 | Li32 convention vs OCT5D2 (their 1/x bug) | NUMERICAL | 3.34e-423 | 423 | atoms step; verify D |
| 15 | duplication Li5(b^4) | NUMERICAL | 2.81e-1629 | 1629 | atoms step; verify D |
| 16 | L(chi8,5) vs M4e | NUMERICAL | 7.09e-1652 | 1652 | atoms step; verify D |
| 17 | local law C_4^(5), C_4^(2,3) m=4 | exact (QS2) | values printed | 50 | local-law step; verify E |
| 18 | verify_oct5d3.py end-to-end | 25 PASS, rc 0 | residuals 1e-349..1e-363 | 360 | verify run, 52.6 s |
"""
open("result.md","w").write(md)
print("result.md complete, %d chars" % len(md))
