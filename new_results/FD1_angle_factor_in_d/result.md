# FD1: the corner-angle factor F_d(c) as a function of d

**Workdir:** `new_results/FD1_angle_factor_in_d/`. Steps under `nice -n 10`, one process at a time; `progress.log` records finished steps.

## Verdict

- **Task B (more data): DONE.** `F_9` and `F_10` are now **EXACT** at all 8 rational angles (single-valued corner/squircle ratio at every cell, polynomial fit in `t = tan^2(a/2)` verified on 9 points each with 3-4 spares). This extends ALPHA2's `d <= 8` to **`d <= 10`**.
- **Side result (squircle):** LEAD1's zigzag top-table closed form `D_d = (-1)^{j+1} binom(d,j) delta_d` is validated at **`d = 10, 11`** (46 cells, 0 mismatches vs a fresh QSI-field squircle computation to `k = 22`); the apparent `d = 10` break was a swapped `A_10/A_11` in ALPHA2's notes, now corrected (`A_10 = 50521`, `A_11 = 353792`).
- **Task A (pattern in d): PARTIAL — CONJECTURE.** No closed form / generating function / low-order recurrence found despite extensive search (Riccati, linear ODE, ratio-of-trig, Hoffman-type families all *excluded*). Several **exact sub-patterns** are established and extended to `d <= 10` (below). One ALPHA2 sub-pattern (`g_3(d) = (14d-81)/315`) is **falsified at d = 9**.
- **Task C (proof): NOT completed.** The corner top recursion (TOW2 port) was analyzed but not implemented in the available time; the value-equation mechanism that would force `F_d` is identified (below) but the proof is open.

## The F_d data (EXACT, d <= 10; F_d = f_d/f_d(1), t = tan^2(a/2))

```
f_2  = t
f_3  = 1 + 3t
f_4  = 2t + 3t^2
f_5  = 2 + 15t + 15t^2
f_6  = 17t + 60t^2 + 45t^3
f_7  = 17 + 231t + 525t^2 + 315t^3
f_8  = 62t + 378t^2 + 630t^3 + 315t^4
f_9  = 62 + 1320t + 5040t^2 + 6615t^3 + 2835t^4     f_9(1) = 15872 = 2*A_9
f_10 = 1382t + 12720t^2 + 34965t^3 + 37800t^4 + 14175t^5   f_10(1) = 101042 = 2*A_10
```
Equivalent ca-form `F_d(ca) = P_d(ca)/(A_d (1+ca)^floor(d/2))`, `P_d in Z[ca]`, `P_d(0) = A_d`:
```
P_2=1-x; P_3=2-x; P_4=5-6x+x^2; P_5=16-13x+x^2; P_6=61-89x+29x^2-x^3;
P_7=272-297x+60x^2-x^3; P_8=1385-2396x+1134x^2-124x^3+x^4;
P_9=7936-10841x+3651x^2-251x^3+x^4; P_10=50521-101187x+62372x^2-12212x^3+507x^4-x^5
```

## Sub-patterns (EXACT on d<=10 data, CONJECTURE in d)

| pattern | status |
|---|---|
| `LC(f_d) = d!/2^{d - s_2(d)}` (s_2 = binary digit sum) | holds d<=10 |
| `f_d(1) = A_d * 2^{s_2(d)-1}` | holds d<=10 |
| `a(2m,1) = a(2m+1,0)` in `f_d = sum a(d,j) t^j` (values 1,2,17,62,1382) | holds d<=10 |
| `P_d(-1) = d!/2^floor((d-1)/2)` | holds d<=10 |
| `P_{2m+1}(1) = A_{2m+1}/2^m`; `P_{2m}(1) = 0` | holds d<=10 |
| `p(2m,m-1) = 2^{2m-1}-m`; `p(2m+1,m-1) = 4^m-m-1` (abs. coeffs of P_d) | holds d<=10 (507, 251) |
| `h_d = f_d/LC`: `[t^{m-1}]h_d = (d-2)/3`; `[t^{m-2}]h_d = (d-4)(5d-13)/90` | holds d<=10 |
| **`[t^{m-3}]h_d`: ALPHA2's `(d-6)(14d-81)/315` FAILS at d=9** (true 88/189, not 3/7) | falsified |

Specializations of the EGF `G(y;x) = sum P_d(x) y^d/d!`: `G(y;0) = sec y + tan y`, `G(y;1) = 1 + sqrt2 tan(y/sqrt2)`, and (if `P_d(-1)` persists) `G(y;-1) = (2+2y+y^2)/(2-y^2)` — these three do not lie in any low-parametric family I tested.

## Method (cheap exact route)

`F_d` needs the top `[k^{d-1}]` of residue polynomials at real poles `{0,+-1,+-r}`, requiring kernels only to `k = 2d`. I extended ALPHA2's exact kernels to `k = 20` (d<=10) and `k = 22` (d<=11) for all 8 angles (~35 CPU-min total, within budget), and the squircle (QSI field) to `k = 22`. Ratios vs the closed-form squircle table are single-valued at every cell. Window-stability of the d=10/11 tops verified across k-windows.

## Files

- `fd1_lib.py` (core: squircle closed-form D-table generator, validated d<=11; zigzag list corrected), `Fd_data.pkl`, `Fd_data_d9.pkl`, `Fd_data_d10.pkl`, `Fd_data_d11.pkl` (F_11 extracted single-valued at 8 angles; final fit interrupted by a file-mode typo in the last step — rerun `solve_poly` on that pkl to print f_11).
- `kern20_*.pkl` (corner, k<=20), `kern22_*.pkl` (corner, k<=22), `kern_sq_11_22.pkl` (squircle QSI, k=11..22), `squircle table generator` in fd1_lib.
- `progress.log`, step files (auto-saved).

**Bottom line:** F_9, F_10 are EXACT (verdict for Task B: done, budget met). The global d-pattern remains **CONJECTURE with d<=10 exact support**; the zigzag numerator of the squircle top table is intact through d=11.