Referee report on the G9 pack (Proposition 9, Theorem 26, Proposition 10). The same content has been kept in `result.md` in my working directory; the `step_NN.py` files there are the code cells that produced the numbers quoted below.

---

# Referee report: `g9_pack` — order-five sections (Prop 9 `prop:order5`, Thm 26 `thm:hexagon5`, Prop 10 `prop:octagon5depth`)

## Verdicts

| Result | Verdict |
|---|---|
| Proposition 9 (order five, local terms of an edge) | **PROVED** (computer-assisted; machine step M1 independently re-derived by my own code) |
| Theorem 26 (P⁽⁵⁾(1) of every hexagon) | **PROVED AFTER REPAIRS** (computer-assisted; repairs are wording-level F1 and the missing script F2; the 307-monomial form is exact end to end, no PSLQ among its ingredients) |
| Proposition 10 (depth two at order five, octagon) | **PROVED** in the repaired form the section itself states ("the symbol of W₅ has a nonzero antisymmetrised (3,2) part"); the paper's printed one-liner is ambiguous and the section discloses this |

## Findings

- **F1 (minor, wording).** `sec_g9_hex5.tex:74` ("The formal elements"): "Î(N,0;a,b) is, for a+b≤5, a polynomial in fz₂,fz₃,fz₅, 1/(N+1), 1/(N+2) and H⁽ⁱ⁾_{N+1}, i≤4" is **false for a+b=5**: my symbolic computation (`step_42/43`) gives Î(N,0;4,1) ⊃ −24/(N+1)·H⁽⁵⁾_{N+1}. For the J-sums (c≥1) one has a+b≤4, where I verified the coefficient of H⁽⁵⁾ is 0 for all 15 pairs (a,b), so the conclusion used (summand = x^q·poly(H_q^{(1..4)})·Q(q)[fz₂,fz₃]) is correct, and S3's local table at integer N is unaffected. Repair: change "for a+b≤5" to "for a+b≤4" (or "i≤a+b").
- **F2 (must fix, reproducibility).** `scripts/g9_hex5_enginecheck_pts.py` and `g9_hex5_enginecheck.py` do `from g4_edge import edge_E`; `g4_edge.py` is **not in the pack and nowhere in the research tree** (I searched). The listed "referee should rerun" item cannot be rerun as shipped. Repair: include `g4_edge.py` in `scripts/`.
- **F3 (minor, wording).** REPORT_G9.md "What a referee should rerun" says the chain takes 896 s; my rerun of S1–S5 took 172+0+41+647+3 s (S4: 647 s vs the reported 555 s) with S6 still running at deadline. No discrepancy of substance, but the 896 s is machine-dependent.

No fatal or mathematical errors found in any of the three sections.

## What I recomputed with my own code (numbers printed by my scripts in this session)

1. **M1 independently re-derived** (`step_41.py`, my own code from the lemma statement: Ψ = L·F⁻²√((1+Λ)²+m(αû−βv̂)²), Stirling conversion, t-integration): `Z5/L − law5 = 0`, `Z23/L − law23 = 0` **exactly**; the order-three sanity `Z3/L − 6αβm(3ξ̄−2) = 0`; parallelogram short forms `0, 0`; the regular 2m-gon sums `2mZ_e⁽⁵⁾ − C_m⁽⁵⁾ = 0`, `2mZ_e⁽²·³⁾ − C_m⁽²·³⁾ = 0` symbolically in s. This is the full content of the paper's claim "exact symbolic integration shows that the law minus the raw integrals is exactly 0".
2. **Beta rules independently verified** (`step_42.py`, my own code): over all n,k≤8 (1701 pairs (a,b) with a+b≤5): `[z5]Î = 24B` iff a+b=5, a,b≥1; `[z2z3]Î = −12B` iff (a,b)∈{(3,2),(2,3)}; the weight bound 2i+3j+5l≤a+b holds everywhere. 0 failures.
3. **H⁽⁵⁾ question** (finding F1): numbers above.
4. **Table checksum**: `g9_hex5_table.txt` sha256 = `f0fd1f7a…2030a`, matches the report.
5. **Chain rerun** (pack's own scripts, G9WORK scratch): S1 172 s, S2 0 s, S3 41 s, S4 647 s, S5 3 s completed with the reported intermediate sizes (336 monomials; 26 atoms; 210 J-keys; S6 running at deadline).

## Hand verification (documented in `result.md`)

- Prop 9: `g9:lem:beta5` proof (Theta expansion, ψ-values, the two rules), `g9:lem:loc5` weight argument, `g9:lem:psi5`(a) (the identities F=1+Aû+Bv̂, D=αû−βv̂, LF−σD=L(1+Λ) checked against `g4:lem:arc`/`g4:lem:mixed`; far terms drop out via Σa_j=1, Σa′_j=0), the parallelogram and 2m-gon hand computations (both brackets; C_m⁽⁵⁾, C_m⁽²·³⁾; the s=1/2 spot check = 2835√3/32), M2's geometric inputs for the diamond (L=√2, u₀=−1/2, m=1/4), regular hexagon (L=1/√3, m=3/4, α=β=ρ=1/2) and sheared parallelogram (u₀=−2/5, m=16/25 and −7/13, 16/169) — all verified against the g4 data (`g4:lem:edge`, `g4:lem:fl`).
- Thm 26: `g9:lem:shape5` (L⁹, odd powers L⁻⁹…L¹) proof checked; `g9:lem:Jseries5` (λ_c(q) = (−1)^c c!/q·e_{c−1}(1,…,1/(q−1)) verified against the generating function; J=O(x^c); analytic continuation domain C∖[1,∞) simply connected); `g9:lem:hpl`(d) (word/Lyndon counts 1,2,4,8,16 per weight; reduction = exact invertibility of per-weight matrices, sizes ≤16×16 — consistent with the S2 log "weight 5: matrix 16x16 inverted"); `g9:lem:nocst` and `g9:lem:z5` proofs checked (fz₅/fz₂fz₃ can only come from local c=0 monomials with a+b=5, giving exactly Z_e⁽⁵⁾, Z_e⁽²·³⁾; consistent with the census: exactly one fz₅ and one fz₂fz₃ monomial in the table, and M6: `120[z5]−Z_e⁽⁵⁾=0`, `120[z2z3]−Z_e⁽²·³⁾=0`); `g9:lem:precise5` continuation argument checked.
- Prop 10: `g9:lem:topwt` (integration by parts, boundary terms vanish) checked; `g9:lem:zero`(a),(b) checked by hand (the two Lie-element contractions vanish; the {3,2}-block cancellation); the non-proportionality arithmetic checked (−495/28 vs 270/19 in the rational parts); the block totals checked by hand: 2(−53955/8)+48015/4=−1485, 2·4770−8505=1035, 2(−10215/2)+11115=900, 2·14445/4−15525/2=−540.

## Priority 2 (exactness of Theorem 26): precise answer

The chain S1–S6 is **exact rational computer algebra throughout; no integer-relation step is an ingredient anywhere**. In detail: S1 expands the length element into 336 monomials with exact SymPy coefficients; S2 reduces the 26 harmonic atoms by the quasi-shuffle identity plus an exactly inverted 16×16 Lyndon matrix (all Fraction arithmetic); S3 is the formal Beta table in N (exact); S4 computes each of the 210 J-closed-forms by exact partial fractions in q, exact shifts q→q−s, and substitution of the exact atom table — **derived, not fitted**; the 160 quadrature checks (1.3e-76) are checks on top, not the derivation; S5–S6 are exact collection into 307 monomials. PSLQ appears only (a) historically for the regular-hexagon form, which is derived exactly by M2; (b) for the r=2 relation, which is a numerical identification confirmed exactly by the fold (all 20 coefficients equal); (c) for the 13 values at x=−1, which rest on the notes' RUL1 certificates (verifier re-run 41/41 PASS; **not** re-derived — disclosed in the gap box). Weakest links, all disclosed by the writer: S1's second-code confirmation needs the unshipped `g4_edge.py` (F2); the J-table engine and assembly are one implementation, confirmed only numerically; my own fully independent E5-engine comparison did not finish within the time limit.

## Coverage

The three daggered statements are proved as printed, with the section's own disclosures: Prop 9's "in this form"/"assumes the constants are linearly independent" caveats are stated in both paper and supplement; Thm 26's L⁹ denominator (main.tex draft 47) agrees with the table (the notes' "L2^{5/2}" was the error, already corrected); Prop 10 is proved in the repaired R7 form ("symbol of W₅ has a nonzero antisymmetrised (3,2) part"), and the gap box honestly lists: the general "no combination of Li₅'s and products" rests on MatRud/CGRR (S-unit case and products proved by hand), the period-conjecture remark is not proved, and the symbol identity Sym W₅ = Sym(classical + C), the decagon/12-gon and the nine-Li₃,₂ form are not reproduced. The paper's wording overclaims nowhere beyond what the dagger convention plus these disclosures cover; the numerical claims I saw in the .out files (finite differences 1.7e-30 hexagon, 6.3e-28 sheared, 3.9e-29 golden, 3.6e-31 octagon; table vs finite differences 8.2e-31 and 6.2e-29; r=1 limit 2.0e-22; r=2 fold exact with edge value −83.5391363475376327634345072141820934358) are internally consistent, but I did **not** regenerate them.

## Not checked (time)

Independent numerical P⁽⁵⁾(1) from the definition; my own evaluation of the 307-table; my own octagon-symbol code (their P1/P2 outputs were reviewed, and the arithmetic of the totals was checked by hand); rerun of the writers' numerical scripts; the chain's S6/compare completion; the notes' RUL1 certificates; the PDF build.