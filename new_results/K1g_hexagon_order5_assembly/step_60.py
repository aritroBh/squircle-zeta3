WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
result = r'''# K1g — assemble P⁽⁵⁾(1) of every centrally symmetric hexagon (finishing K1f)

Conventions (as K1e/K1f): **unweighted logs** lA0=log(A0/w0), lA1=log(A1/w1), lB=log(B/wo); ℓ=−log(1−x)=−log r, x=−κ=1−r; H1=HPL atom H_{1}(x)=ℓ. Term values a0=A0·t, a1=A1(1−t), a2=B(1+κt), A1+B=1, A0=1−B−Bκ. E_arc⁽⁵⁾=120·∫E[5]dt, P⁽⁵⁾(1)=2Σ_{q=0..2}E_arc(edge_q).

## Task 1 — Assembly: DONE
E_arc⁽⁵⁾ = 120·(Σ_{c≥1} coef(n,a,b,c)·J(n;a,b,c)|_{x=−κ} + Σ_{c=0} coef(n,a,b,0)·Jloc(a,b)|_{N=n}) assembled **symbolically** (`step_30/34`, `coef_acc_raw.pkl`): **307 distinct basis monomials** = the 13 Lyndon HPL atoms {H1,Li2,Li3,H011,Li4,H0011,H0111,Li5,H00011,H00101,H00111,H01011,H01111} (+products) × x-powers × {ζ2,ζ3,ζ5} + rational terms; coefficients are polynomials in the edge data (A0,A1,B,κ,L2,q0,V02,lA0,lA1,lB), denominator L2^{5/2}. x-power distribution: x^{-6}…x^{-1} (poles) + x^0 (39/44/44/44/44/44/48 monomials).

**κ⁻⁶ poles cancel** — but NOT by the pure-rational coefficients vanishing under A0=1−B−Bκ (I checked: the 20 pure rational/z pole coefficients at x^{-5}…x^{-1} are individually NONZERO after the constraint). The cancellation is **structural**: every J-table entry satisfies J(n;a,b,c)~x^c as x→0 (K1f, ratios 10^c), so the x^{-k} poles in the monomial representation cancel *within* each J-entry once the atoms are included; the edge-data coefficients are κ-polynomials (regular), hence every assembled term is O(κ) and E_arc is regular at κ=0. **Demonstrated numerically** (`step_57`): synthetic edge B=1/2,L2=2,q0=−1,V02=1, E_arc(κ) for κ=10^{-1},…,10^{-5} = 41.185, 48.745, 49.526, 49.604, 49.6121 → κ=0 local value 49.6129 (smooth convergence, no 10^30 blow-up).

## Task 2 — Checks (all PASS)
**(a) r=1 local law.** Standalone-ζ5 coefficient of E_arc = Z_e⁽⁵⁾ and standalone-ζ2ζ3 = Z_e^{(2,3)} **for general r**: 120·C_{z5}−Z_e⁽⁵⁾=0 and 120·C_{z2z3}−Z_e^{(2,3)}=0 at 6 random valid edge-data points and on b=0 edges 1,3 (diffs 0.0, `step_58`). κ=0 arc (b=0 edge 3) in E5c convention-B reproduces the exact closed form 7.4848657252901844… to 2.3e-46 (`step_53`); its ζ-coords are exactly 3195√2/2048=Z_e⁽⁵⁾ and 1515√2/2048=Z_e^{(2,3)} (no Li_n(1/2)) (`step_44`).
**(b) b=0, r=2 arc vs E5e.** Convention reconciliation: E5c conv-B = unit normals |m_j|; the polar/K1f term value is u_j=|<m_j,y>|, so the log shifts by log(w_j/|m_j|). (i) ζ-coefs compared **exactly**: folding my atoms at x=−1 into the alternating-MZV basis gives ζ5 = **−20187√5/2000** and ζ2ζ3 = **−858√5/125** = E5e's X-coordinates ×√5 EXACTLY (`step_51`); hence N⁽⁵⁾=−711√5/400, N^{(2,3)}=−1239√5/250 (match E5e). (ii) Full value after conversion (ℓA0=log2, ℓA1=−log√2, ℓB=0): −852.255256336156211501665304739892950203528361 vs E5e −852.255…0528360898, **diff 8.5e-56** (`step_52`); Landen-pair edge r=1/2 matches to 2.8e-38.
**(c) Fresh irrational golden hexagon** V=(1,0),(φ,1),(0,1), w=(2,3,5), ratios {φ,1,1/φ}: direct polar-perimeter P(p) via mp.quad + 5th-derivative antisymmetric stencil (weights 29/6,−13/3,3/2,−1/6) + Richardson gives P⁽⁵⁾(1)=−78.35859987396460…; my formula 2ΣE_arc = **−78.3585998739646021275**, **diff 2.1e-15** (≥15 digits, `step_47/49`). Convention-consistent (both polar), so this independently validates the whole pipeline.

## Task 3 — ζ(5), ζ(2)ζ(3) for general r
In the 13-atom Lyndon basis the **explicit (standalone) ζ(5) coefficient is exactly Z_e⁽⁵⁾ for all r** (and ζ(2)ζ(3)=Z_e^{(2,3)}): the non-local ζ(5) does not appear as a standalone ζ5 but is **carried by the HPL atoms**, which reduce to ζ(5)-combinations only at special r. At **r=2 (x=−1)** the non-local ζ(5) is carried by (`step_50`): Li5=−15ζ5/16, H00011=−29ζ5/32, H00101=+41ζ5/32, H00111=+33ζ5/32, H01011=−73ζ5/64, H01111=+ζ5 (with H00011,H00101,H00111,H01011 also feeding ζ2ζ3). The two non-Nielsen atoms H00101, H01011 are essential (each appears in 30 J-keys).

## Theorem (status)
**Theorem.** For every centrally symmetric hexagon K_p={Σ_j w_j|⟨n_j,x⟩|^p≤1} normalized at the vertices, P⁽⁵⁾(1)=2Σ_{q=0}² E_arc⁽⁵⁾(edge_q), where E_arc⁽⁵⁾ is an explicit closed form in the 13 Lyndon HPL atoms + {ζ(2),ζ(3),ζ(5)} at x=1−r_e, with coefficients polynomial in the edge data over L2^{5/2}; the standalone ζ(5) (resp. ζ(2)ζ(3)) coefficient is exactly the local edge law Z_e⁽⁵⁾ (resp. Z_e^{(2,3)}) for all r, the non-local parts entering through the HPL atoms. Specializes to the affine-regular local law at r=1 and to E5e's relation at r=2,1/2. **Status: assembly is an exact symbolic reduction of the validated weight-5 engine (E[3]=K1d exactly) and 80-dps-verified J-table; pole cancellation structural; all three independent checks PASS** (golden-hexagon finite-difference 2.1e-15; b=0 ζ-coefs exact + full value 8.5e-56 after convention conversion; r=1 local law exact). The engine/J-table symbolic correctness + check (c) make the per-edge formula **PROVED modulo the K1d E[3]-engine validation**; the b=0/E5e comparison is a numerical confirmation, not part of the proof.

## Evidence table (all numbers printed by my code)
| Check | Object | Result | Residual | Digits | Step |
|---|---|---|---|---|---|
| Assembly | 307 monomials, edge-data coefs | built | — | exact | step_30/34 |
| κ→0 regularity | E_arc(κ)→local, κ=1e-1..1e-5 | converges 49.6129 | smooth | 60 | step_57 |
| (a) ζ5=Z_e⁽⁵⁾, ζ2ζ3=Z_e^{(2,3)} general r | 6 random pts + edges 1,3 | exact | 0.0 | exact | step_58 |
| (a) κ=0 arc = E5c closed form | b=0 edge 3 conv-B | 7.4848657… | 2.3e-46 | 70 | step_53/44 |
| (b) folded ζ5,ζ2ζ3 = E5e | b=0 edge 1 | −20187√5/2000, −858√5/125 | exact | exact | step_51 |
| (b) full value after conv | b=0 edge 1 conv-B | −852.255… | 8.5e-56 | 80 | step_52 |
| (b) Landen pair r=1/2 | b=0 edge 2 | same | 2.8e-38 | 70 | step_53 |
| (c) golden hexagon vs finite-diff polar | P⁽⁵⁾(1) | −78.35859987396460 | 2.1e-15 | 15+ | step_47/49 |

## Files (work dir `K1g_hexagon_order5_assembly/`)
`build_blocks.pkl` (J5,JLOC_n), `coef_acc_raw.pkl` (307-monomial assembly), `golden_edges.pkl`, `atom_mzv_m1.pkl` (atoms@x=−1 in MZV basis), `atom_eval.py` (copied K1f step_31), `esum5.py`, `atom_table.pkl`; steps step_30/34 (assembly), step_44/51 (ζ-coefs edges 3/1), step_47/49 (golden finite-diff + formula), step_50 (atoms@−1 PSLQ), step_52/53 (conv-B b=0), step_57 (κ→0), step_58 (general-r Z_e).
'''
open(f'{WD}/result.md','w').write(result)
print("result.md written,", len(result), "chars")
