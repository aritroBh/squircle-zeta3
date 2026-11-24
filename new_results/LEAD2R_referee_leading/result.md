# LEAD2R — referee report on LEAD2's leading-coefficient proof

**Scope.** Re-derived every step of `proof_lead.md` by hand and re-checked every identity with my own code (`lead2r_lib.py`, written from the proof text — no import of `lead2_lib`), at 50–110 digits. Verdicts below; residuals are all printed numbers from my own runs.

## 1. Verdict table (Task 1 items)

| item | verdict | evidence (my code, my residuals) |
|---|---|---|
| Lemma 2.1 (u-substitution, 4 measures, `m_odd=1/(1-u)-u m_e`, parity, endpoint limits −1/4, 1/4, deck `v_+v_-=−1`, pole pairing) | **CORRECT** | sympy `radsimp==0` on all 16 identities, incl. `m_oddlast~1/(2(1-u))` at u=1 and finiteness at u=−1 |
| Eq. (3) + convergence of all M-integrals | **CORRECT, one documentation typo** | full key sum with 4 measures (incl. diagonal/oddlast) `= (d!/2)e_{d+1}` for d=1..13 both parities, worst 1.1e-24 (dps 60, precision-limited trend; d=13 assembly 1.1e-65 at dps 110). **The printed j-range `0..floor((d-1)/2)` omits the even-d diagonal key** — using it alone gives errors 0.26–1.8 (printed); Lemma 4.1's proof adds it back correctly. MINOR FIX: change the range to `0..floor(d/2)` |
| Master Lemma (i) WW=E+E−, (ii) FF=(E−E−)/s, (iii) GG=2π/sin−(E+E−)/s; Beta form (4); IBP; analytic continuation; signs | **CORRECT** | re-derived by hand (checked every substitution, boundary terms vanish exactly on −1<Re s<1; the two-piece split `GG(s)=M_g(−s)−M_g(1−s)` is convergent on the whole strip). Numerics: 10 strip points incl. s=±0.99, ±0.9, 0.2+0.4i, −0.35−0.28i, 0.01+0.9i, −0.97+0.15i at dps 55: WW worst 3.6e-41; FF/GG/IBP worst 7.4e-26, with the F partial-fraction series agreeing with quadrature to the same 7.4e-26 (two independent methods). Inversion symmetries `h(1/w)=w²h(w)`, `q(1/w)=w²q(w)`, `g1(1/w)=−w²g1(w)`: symbolic True |
| F-series erratum to LEAD1 | **CORRECT** | LEAD1's printed form (both sums from n≥1) is off by exactly `−1/(1−y)`: measured −1.43 vs −1/(1−y)=−1.43 at y=0.3; −0.69 at y=−0.45. Corrected form (second sum from n=0) matches quadrature 3.1e-39 / 6.7e-52 |
| Lemma 5.1 (`u m_e` collapse, GF=GG/4, eq. (7)) | **CORRECT** | parity/collapse re-derived; (7) numerically for odd d=1..13: worst 1.9e-15 (dps 55) → 1.1e-65 (dps 110, d=13); growth was cancellation, not error |
| Splitting `m_odd=1/(1-u)−u m_e` termwise legitimate | **CORRECT** | every piece separately convergent (each key has i=d−j≥1, so `L1^i/(1-u)` is O((1-u)^{i-1}) at u=1; `u m_e` continuous on [−1,1]); verified by evaluating each piece's integrals separately |
| Lemma 5.2 (`H(x,y)=B(1+x,y)`, Laurent bookkeeping) | **CORRECT** | integral identity at 6 points incl. the coordinator's FAIL point (0.3,0.4): 2.1e-41 at dps 55 (the coordinator's FAIL was a quadrature tolerance in their script, not a math issue). Laurent bookkeeping: `aa_{p,q}` from the integral definition == `[x^p y^{q+1}]P`, p=1..6,q=0..5, worst 8.9e-34 (dps 55); `P(x,0)=P(0,y)=1` exact ⇒ `aa_{0,d}=0` for d≥1, legitimately extending the sum to p=0. Note: the displayed series equality is best read coefficientwise (the proof's own analytic argument — `a_p(y)` analytic on Re y>−1 for p≥1 — is the valid content); the double series converges for \|x\|<1 |
| Lemma 5.3 (`yB(1+x,y)=xB(1+y,x)`, recursion (10), ranges p≥1,q≥0) | **CORRECT** | recursion `aa_{p,q}=aa_{q+1,p-1}` verified EXACTLY (symbolic series of P in γ_E, ζ(n)) for p≤8, q≤7 |
| Theorem 5.4 (involution σ(p,q)=(q+1,p−1) on S_d, fixed point D_*, σ(A_*)=(S_d\A_*)∪{D_*}, halved last key, sign of B(1+z,−z)=−π/sin, value (1−2^{−d})ζ(d+1)) | **CORRECT** | fold identity `sum_{A_*}f−½f(D_*)=½sum_{S_d}f` verified EXACTLY for odd d=1..13; values printed: π²/12, 7π⁴/720, 31π⁶/30240, 127π⁸/1209600, 73π¹⁰/6842880, 1414477π¹²/1307674368000, 8191π¹⁴/74724249600 = (1−2^{−d})ζ(d+1). Sign check: `B(1+z,−z)=Γ(1+z)Γ(−z)=−π/sin(πz)` confirmed (also numerically, residual 0.0). T_d via the defining bracket (9) by quadrature: d=1..13, worst 3.7e-15 (dps 55) → 2.1e-65 (dps 110) |
| Theorem 5.5 (exact cancellation of the two (1−2^{−d})ζ(d+1) finite parts) | **CORRECT** | assembly `T_d − (u m_e part) = (d!/2)e_{d+1}`: odd d=1..13, worst 1.9e-15 → 1.1e-65 (dps 110) |
| Lemma 4.1 / Theorem 4.2 (even-d parity collapse, single integral) | **CORRECT** | collapse worst 5.6e-54, middle-key correction `M(d/2,d/2;m_diag−m_e/2)=0` by oddness (numerically ≤1.6e-63), single-integral form worst 2.9e-22 (d=2..12) |
| Lemma 6.1 / Theorem 6.2 / Corollary 6.3 (zigzag signs, assembly, Hadamard GF, duplication, log E) | **CORRECT** | `δ_d=eps_{d−1}A_d/d!` and `eps_mA_{m+1}/m!=[z^m](1−sinh z)/cosh²z` symbolic for d,m≤21 (two independent zigzag computations agree); E-coefficients by log E vs direct Gamma-product Taylor agree to 2.3e-61 |

## 2. What rests on the table (Task 2) — with one status UPGRADE

I re-extracted the top coefficients myself from the raw kernels (`WTS kern_1..18.pkl` + `TWOEND1/red/kern_19..21.pkl`), exact rational arithmetic:

- **(a) Polynomial law K(d)** (kernel residues polynomial in k of degree d−1 for k≥d+1): exact data (d-th differences vanish identically) for **d=1..10**, spare Newton points: d≤8 spare ≥5, **d=9 spare = 3** (LEAD2 had ZERO spare using k≤18; TWOEND1's k=19,20,21 upgrade it), **d=10 spare = 1** (LEAD2 listed d≥10 as hypothesis). For all d: OPEN (= tower law T; TOW2 has no result yet). TOWN1's exact towers (depths ≤8, levels ≤52) + WTS Lemma 3 cover d≤5 to level 52. Note: proof_lead.md §1 says "polynomial of degree d / d-th differences constant" — the code and the rest of the document use degree d−1 ((d−1)-st differences constant = Dtop, d-th vanish); wording slip, MINOR FIX.
- **(b) D-table values** (binomial × zigzag): matches my extraction at **every** key and pole for d=1..10, including the halved odd-d last key (pole a) and the odd-d diagonal carrying −b and 0 halved. Exact data d≤10; hypothesis d≥11.
- **(c) Structural hypotheses**: all √2-parts of top coefficients vanish (d=1..10, exact); pole set is exactly {0, ±a, −b} (+b never occurs, no ±i poles, no constant part) — exact data k≤21.

**What is proved for which w:** given the D-table at d=w−1, the proof establishes the law for **every** w (Task-1 chain, all verified). As data the table is now exact for d≤10 ⇒ law holds for the polynomials fitted to exact level data for **w≤11** (w=11 = d=10 with 1 spare point). "Proved for w≤10" is defensible but understates/overstates: with spare-point exact data the range is w≤11; strictly unconditional (all-k polynomiality proved, not finite data) only w=2 (d=1) via WT2b/CX5, plus w=3 via TOW1's proved depth-≤2 tower law + Lemma 3; w≤6 has exact tower-data support to level 52 via TOWN1.

## 3. End-to-end (Task 4)

My own pipeline (own D-extraction incl. the new d=10 row; own N-integral quadrature at dps 80): printed `|c(w)−law| = |law|` for every w=2..11 (0.531, 0.463, 0.512, 0.409, 0.336, 0.253, 0.189, 0.137, 0.0985, 0.0696). This pattern is exactly `|2·law − law|`: my slice loop integrated [b,a] twice (verified at d=1: single-loop N(−a)=−0.374929952811153 equals plain quad and numpy to 1e-6+). So my computed sum equals **2·c(w) with c(w)=law to printed precision** for w=2..11 — consistent with the law at ≥3 digits plus every ingredient at 40+ digits; the corrected residual table was not printed before the turn limit. LEAD2's own verify_lead2.py part F (kernels→c(w), w=2..10, residual 8.5e-91) was rerun PASS by the coordinator; all of its inputs (D values, measures, law coefficients, quadrature identities) are independently confirmed here. The end-to-end numeric confirmation is therefore credible but my independent full-precision rerun remains **incomplete** (the only unfinished item).

## 4. Final verdicts (Task 5)

| statement | verdict |
|---|---|
| (i) even-weight step (odd d) given the D-table | **PROVED** (involution, finite-part cancellation, Master Lemma all verified by hand + exactly/numerically) |
| (ii) the law for all w given the D-table for all d | **PROVED** (both parities; LEAD1 even-d part refereed: OK, one docs typo in the printed F-series, no mathematical gap) |
| (iii) unconditional law for small w | **PROVED AFTER REPAIRS** in the sense: law holds for the exact-data polynomials for w≤11 (d≤10, spare points, now incl. d=9,10); strictly all-k unconditional only w=2 (and w=3 via TOW1). Paper wording must keep the condition explicit |

Minor fixes to record: (1) eq. (3) j-range should be `j=0..⌊d/2⌋` with the measure dictionary (diagonal key currently outside the printed range); (2) §1 "polynomial of degree d" → "degree d−1"; (3) Lemma 5.2's series equality should be stated coefficientwise (the proof's analytic argument already is); (4) the erratum on LEAD1's F-series (second sum from n=0) is correct and already propagated to the paper.

**Proposed Theorem wording (paper).** *Let c(w)=[k^{w−2}]W_w(k). Assume the degree-(w−1) kernel residues are polynomials in k of degree ≤w−2 for k≥w with top-coefficient table β_d(r) given by the binomial–zigzag pattern (exact data for d≤10 — verified with spare Newton points through level 21, and through level 52 for d≤5 via TOWN1; implied in general by the tower law, Conjecture X). Then for every w≥2 for which the table holds at d=w−1, c(w)=σ_{w−2}A_{w−1}/(2(w−2)!)·[y^w]𝓔(y), 𝓔(y)=Γ(1+y/2)Γ(1/2−y/2)/√π; equivalently Σc(w)u^w = ½[u²(1−sinh u)/cosh²u ⊙ 𝓔(u)].*

## 5. Files produced (this folder)

- `lead2r_lib.py` — independent library (zigzag two ways, δ_d, E-coefficients two ways, kernel loaders, exact Newton differences)
- `lead2r_Dtable.pkl`, `lead2r_cw_resid.pkl` — extracted D-table d=1..10 with spare counts; end-to-end residuals
- `progress.log`; step files (harness-saved): D-table extraction+structural checks (step_04/05), Lemma 2.1 symbolic (step_05), Master Lemma + erratum (step_06–13), Lemma 5.2/H=B + Laurent bookkeeping (step_14–18), T_d/um/assembly d=1..13 + dps-110 reruns (step_19–21), even-d collapse + eq(3) j-range demonstration (step_22b), end-to-end (step_23–24 debug). A consolidated `verify_lead2r.py` was not assembled within the turn limit; each check above is a saved step file with printed residuals (all PASS except my own factor-2 slice bug in step_23, documented and root-caused above).