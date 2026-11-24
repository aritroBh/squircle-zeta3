result_md = r'''# DOD5D — the regular 12-gon's P⁽⁵⁾(1) needs a depth-two polylogarithm: RESULT

## Bottom line

**YES — like the octagon, the regular 12-gon needs a depth-two polylogarithm at order five**
(PROVED at the symbol level; number level modulo the period conjecture). With
\(u=\sqrt{2-\sqrt3}=(\sqrt6-\sqrt2)/2\) and \(V=\langle X=\log(2+\sqrt3),\;Y=\log(1+\sqrt3),\;Z=\log\sqrt3\rangle\)
(\(\dim V=3\)), the antisymmetrised (3,2) cobracket of the weight-5 part of \(P^{(5)}(1)/u\),
paired against the Hall bases of \(\mathrm{Lie}_3(V)\otimes\mathrm{Lie}_2(V)\) (8×3 = 24 pairings,
M4c's exact convention), is **nonzero in all 24 pairings**. For example
\[
\bigl\langle [X,[X,Y]]\otimes[X,Y]\bigr\rangle \;=\; 2160\,(66-37\sqrt3)\;\neq\;0 .
\]
The obstruction is **exactly** the cobracket of an explicit 24-term depth-two combination
(21 \(\mathrm{Li}_{4,1}\) + 3 \(\mathrm{Li}_{3,2}\) at \(S\)-units of \(\mathbb Q(\sqrt3)\),
residual exactly 0): after subtracting it, the symbol is classical (depth one), i.e.
\(S(W_5/u)=S(\text{classical}+\text{depth-two})\). All four controls pass: \(\mathrm{Li}_5\) at 13
field \(S\)-units and 300 random shuffle products pair to 0; M4c's \(W_5\) reproduces
\((84,-38\sqrt2)\); OCT5D2's octagon reproduces \((-1485+1035\sqrt2,\;900-540\sqrt2)\) with the
same general code.

## Task 1 — Value (NUMERICAL, two quadrature levels)

One-wedge ×12 engine (`dod_engine.py`, adapted from OCT5D's validated oct_engine): tanh-sinh,
dps 120 (level diff \(9.10\times10^{-121}\)) and dps 170 (level diff \(7.53\times10^{-172}\)); the two runs
agree to \(1.31\times10^{-117}\).
\[
P^{(5)}(1)=-62.53724775498639535309129282205803572684769277410373262688145818074695731422874528214220721559218128482460652943087962941126065293896549327018536685224048730049659698854\ldots
\]
\[
P^{(5)}(1)/u=-120.812685423158829233550820329306095036347192964810348085935676572312012620107516072347383028089431895314864448521566618\ldots
\]
Sanity: \(P(1)=12\sqrt{26-15\sqrt3}\) (residual \(1.48\times10^{-181}\) at dps 190);
\(P'''(1)\) vs `DOD1_dodecagon/P3_630.pkl` (residual \(2.63\times10^{-174}\)).

## Task 2 — Exact structure (PROVED algebra; NUMERICAL check 140 digits)

Edge (wedge \((0,\pi/6)\), signed normals \(s=(+,+,+,+,-,-)\)): from `symeng.edge_params`:
\(j_0=3,\;j_1=4\), non-local \(=\{0,1,2,5\}\); \(\alpha=\beta=(2-\sqrt3)/2\),
\(L=u(2-\sqrt3)=\sqrt{26-15\sqrt3}\), \(h=u/2\), \(s_0=-L/2\). Four non-local ratios
\[
\kappa=\Bigl(\tfrac{\sqrt3-2}{2},\;\tfrac{2\sqrt3-3}{3},\;\sqrt3-1,\;\tfrac{\sqrt3-3}{3}\Bigr).
\]
Exact substituted \(E_5(t)\) (from `E5_edge_law/step_02_full.pkl`, 1910 terms): every term has
\(L{+}h{+}s_0\) degree exactly 1, so every monomial coefficient is \(u\cdot\mathbb Q(\sqrt3)\)
(checked: no \(\{1,\sqrt3\}\)-part in any of the 11006 monomials — the analogue of the octagon's
\(u=\sqrt{2+\sqrt2}\) factor). The exact expansion matches `symeng.E5_sym` to \(3.3\times10^{-114}\)
(dps 120). 11006 monomials \(t^n(\log t)^a(\log(1-t))^b\prod_i(\log(1+\kappa_i t))^{c_i}\times\)
constant logs (weight distribution 17/134/659/2474/7722 for weights 1..5).
\(P^{(5)}(1)=12\cdot5!\int_0^1E_5\,dt=1440\int E_5\) (quadrature check \(1.5\times10^{-50}\) at dps 60).
Assembled from 2772 numerically integrated monomials (`itab140.pkl`, flint arb tanh-sinh,
level diff \(2.75\times10^{-148}\)):
\[
\Bigl|1440\!\!\sum_{\text{mons}}\!\!\mathrm{coef}\cdot I_{\mathrm{mon}}\;-\;P^{(5)}(1)/u\Bigr|
= 3.55\times10^{-140}\quad(\ge 60\text{ digits required } \checkmark).
\]

## Task 3 — Letters and V (PROVED)

\(G\)-letters: \(\{0,1,p_0,p_1,p_2,p_3\}\) with \(p_i=-1/\kappa_i\):
\(p_0=4+2\sqrt3=V^2,\;p_1=-(3+2\sqrt3)=-\sqrt3\,U,\;p_2=-(1+\sqrt3)/2=-U/V,\;p_3=(3+\sqrt3)/2=\sqrt3\,U/V\)
(\(U=2+\sqrt3,\;V=1+\sqrt3\)). Closed under \(t\mapsto1-t\): \(p_0\leftrightarrow p_1,\;p_2\leftrightarrow p_3\).
**All 30 pairwise letter differences factor in \(\langle U,V,\sqrt3\rangle\)** (exact `logvec`
factorisation, all asserted): e.g. \(p_0-p_1=U^2\), \(p_2-p_3=-U\), \(p_1-p_3=-\sqrt3 UV/2\).
The prime 13 does **not** enter; no other primes appear. Hence
\(V=(S\text{-units}\bmod\text{torsion})\otimes\mathbb Q=\langle X,Y,Z\rangle\), \(\dim V=3\)
(the constant logs are \(-2Y\), \(-X\), \(Z-2Y\); \(\log2=2Y-X\)).
Exact rewrite of all 2772 monomial integrals into shuffle-regularised \(G(w;1)\)
(7776 words = exactly the words over 6 letters not starting with 1 — the octagon's
regularisation consistency): numerical evaluation via the six-letter \(G\)-table
(`gtab6_d130_N320.pkl`, 9115 values; self-checks: quadrature \(2.1\times10^{-132}\),
shuffle \(1.4\times10^{-132}\), log identity \(1.1\times10^{-131}\), polylog identities
\(\le6.5\times10^{-99}\)) assembles \(P^{(5)}(1)/u\) to \(1.13\times10^{-96}\) of the engine.
Weight-5 part of \(P^{(5)}(1)/u\) (identical from itab and gtab assemblies):
\[
-175.610996035664253403647607710697722922183103424161230954083\ldots
= \underbrace{-186.70176165766197128755378942475431792534473368448}_{\text{local}}
+\underbrace{3.5485647202796266786959620545770716523377118289154}_{\text{single}}
\]
\[
\underbrace{+7.4037148383166814368778080221642790464886879597396}_{\text{double}}
+\underbrace{0.12740092441561833764559999527385630891109675603075}_{\text{triple}}
+\underbrace{0.011085138985791430686811642041387995424133715636783}_{\text{quad}} .
\]

## Task 4 — Depth test (PROVED at symbol level)

Full number-level Goncharov symbol with every constant entry \((a_i-a_j)\) (constant-retaining
`g_symbol` over the \(\mathbb Q(\sqrt3)\) \(S\)-unit logvec, same formula as the validated
`A3_W4_proof/w4_symbol.py`). Weight-5 symbol \(S(W_5/u)\): 243 entries (all words over
\(\{0,1,2\}\) of length 5); assembled twice by independent code paths (by-const-pattern grouping
vs per-monomial), identical symbols and obstructions. The antisymmetrised (3,2) tensor
\(\Psi=\sum_w s(w)\,(w_{[:3]}\otimes w_{[3:]}-w_{[2:]}\otimes w_{[:2]})\) is paired against the
Hall bases of \(\mathrm{Lie}_3(V)\otimes\mathrm{Lie}_2(V)\) (this annihilates exactly the shuffle
products by Ree's theorem; M4c's convention, calibrated below).

**Total obstruction table** (entries \(p+q\sqrt3\); \(r3=\sqrt3\); rows \(\mathrm{Lie}_3\), columns \(\mathrm{Lie}_2\)):

| Lie₃ \ Lie₂ | [X,Y] | [X,Z] | [Y,Z] |
|---|---|---|---|
| [X,[X,Y]] | 2160(66−37r3) | 135(209−137r3) | 270(−603+358r3) |
| [Y,[X,Y]] | 540(−357+205r3) | 270(106−55r3) | 2160(27−17r3) |
| [Z,[X,Y]] | 270(119−71r3) | 270(−153+86r3) | 270(235−127r3) |
| [X,[X,Z]] | 135(−733+403r3) | 135(−8211+4730r3) | 4185(575−331r3) |
| [Y,[X,Z]] | 1080(158−89r3) | 135(16179−9311r3) | 270(−17473+10061r3) |
| [Z,[X,Z]] | 135(185−119r3) | 135(−12478+7201r3) | 135(26401−15243r3) |
| [Y,[Y,Z]] | 1080(−274+155r3) | 270(−16566+9533r3) | 540(17807−10254r3) |
| [Z,[Y,Z]] | 270(−171+113r3) | 135(25431−14683r3) | 1080(−6714+3877r3) |

All 24 pairings are nonzero ⇒ **depth two is needed**. Block decomposition (sum = total, checked):
**local block obstruction = 0 exactly** (consistent with `prop:order5`: local terms are
classical); every non-local block — the four single-ratio, the six double, the four triple and
the quadruple block — has nonzero obstruction. Mirror equalities are exact at the obstruction
level: single-0 = single-1, single-2 = single-3, double-02 = double-13, double-03 = double-12,
triple-012 = triple-013, triple-023 = triple-123 (the \(\kappa\)-mirror pairs
\(t\mapsto1-t:\;p_0\leftrightarrow p_1,\;p_2\leftrightarrow p_3\)). Unlike the octagon — where
contrary to the planner heuristic even the single-ratio block carried depth two — here every
mixed block carries it, and so do the triple and quadruple blocks (which the octagon did not have).

**Controls (all PASS, run as OCT5D2 did):** \(\mathrm{Li}_5\) at 13 \(S\)-units of \(\mathbb Q(\sqrt3)\)
pairs to 0; 300 random shuffle products pair to 0; M4c's \(W_5\) (its own letters) pairs to
\((84,-38\sqrt2)\); OCT5D2's octagon block symbols pair to \((0,0)\) local,
\((-1485+1035\sqrt2,\;900-540\sqrt2)\) total — with this same general (any \(\dim V\), any exact
coefficient ring) code. G-symbol identities: \(S\,G(0^{k-1},p_0;1)=-S\,\mathrm{Li}_k(1/p_0)\)
(\(k=2,5\)), \(S\,G(p_0)\,\ш\,S\,G(p_2)=S\,G(p_0,p_2)+S\,G(p_2,p_0)\) — all exact.

**Depth-two correction (PROVED at symbol level):** the obstruction equals, exactly, the
cobracket of (residual 0 in all 48 rational equations, flint fmpq rref rank 48):

```
  (2123550-1222560r3) * Li_4,1( (-6+4r3), (1/2+1/6r3) )
  (-1898640+1092690r3) * Li_4,1( (-5+3r3), (1/2+1/2r3) )
  (-1070955+617085r3) * Li_4,1( (-3+2r3), (1/3r3) )
  (-6197715/2+1785510r3) * Li_4,1( (-3+2r3), (1/2+1/6r3) )
  (-2718765/2+1568025/2r3) * Li_4,1( (-3+2r3), (1/2+1/2r3) )
  (-8088525/2+2326455r3) * Li_4,1( (-3/2+3/2r3), (1/3r3) )
  (-3276855/2+942165r3) * Li_4,1( (-1+2/3r3), (1r3) )
  (-1116315+644490r3) * Li_4,1( (-1+1r3), (-1/2+1/2r3) )
  (3868695-2229120r3) * Li_4,1( (-1/3+1/3r3), (1r3) )
  (-769095/2+224775r3) * Li_4,1( (1/3r3), (-1+1r3) )
  (-4944105/2+1423170r3) * Li_4,1( (1/3r3), (-1/2+1/2r3) )
  (3187755/2-1832355/2r3) * Li_4,1( (1/2r3), (-1+1r3) )
  (4523715/2-1307745r3) * Li_4,1( (1/2r3), (1-1/3r3) )
  (-2472795/2+710370r3) * Li_4,1( (1/3), (3/2-1/2r3) )
  (-6692085/2+1930500r3) * Li_4,1( (1/2-1/6r3), (1r3) )
  (-1953855+1130895r3) * Li_4,1( (1/2), (-1+1r3) )
  (1448145/2-419985r3) * Li_4,1( (1/2), (1+1/3r3) )
  (-5156055/2+1485405r3) * Li_4,1( (1/2+1/6r3), (-1+1r3) )
  (-974430+561735r3) * Li_4,1( (1/2+1/6r3), (2-1r3) )
  (-587655/2+162810r3) * Li_4,1( (2/3), (3/2-1/2r3) )
  (-1046115/2+301725r3) * Li_4,1( (1+1/3r3), (2+1r3) )
  (1643220-948915r3) * Li_3,2( (-1/3+1/3r3), (1r3) )
  (1075680-618435r3) * Li_3,2( (1/2-1/6r3), (3/2+1/2r3) )
  (-1116990+641655r3) * Li_3,2( (1/2+1/6r3), (3-1r3) )
```

(all points \(S\)-units of \(\mathbb Q(\sqrt3)\) with \(1-x\), \(1-y\), \(1-xy\) also \(S\)-units,
verified by `logvec`; \(\mathrm{Li}_{a,b}(x,y)=G(0^{b-1},1/y,0^{a-1},1/(xy);1)\).) The selected
vectors are the 1200 lowest-height of 4487 distinct nonzero cobracket vectors from 1822 valid
pairs × 4 MPL families.

## Evidence table

| # | check | residual / result | digits | producer |
|---|---|---|---|---|
| 1 | \(P^{(5)}(1)\), dps-120 vs dps-170 engines | agree to \(1.31\times10^{-117}\); level diffs \(9.1\times10^{-121}\), \(7.5\times10^{-172}\) | 117/170 | dod_engine.py, step_19/20 |
| 2 | \(P(1)=12\sqrt{26-15\sqrt3}\) | \(1.48\times10^{-181}\) | 181 | step_21 |
| 3 | \(P'''(1)\) vs P3_630.pkl | \(2.63\times10^{-174}\) | 174 | step_21 |
| 4 | edge params (α,β,L,h,s0,j's,κ) exact | printed equalities to 30 digits | 30 | step_23 |
| 5 | exact substituted \(E_5\) vs `E5_sym` | \(3.3\times10^{-114}\) | 114 | dod5d_exact.py, step_36 |
| 6 | \(P^{(5)}=1440\int E_5\) | \(1.46\times10^{-50}\) | 50 | step_37 |
| 7 | 11006 monomials, coeffs all \(u\cdot\mathbb Q(\sqrt3)\) | exact (asserted) | — | dod5d_exact.py, step_35 |
| 8 | monomial integrals itab140 (2772) | level diff \(2.75\times10^{-148}\), ball rad \(1.4\times10^{-147}\) | 148 | itab_build2.py |
| 9 | itab assembly vs engine \(P^{(5)}/u\) | \(3.55\times10^{-140}\) | 140 | dod5d_assembly.py |
| 10 | exact G-rewrite, 2772 monomials, 7776 words | exact; no leading-1 words | — | dod5d_grewrite.py, step_65 |
| 11 | G-table self-checks | \(2.1\times10^{-132}\) (quad), \(1.4\times10^{-132}\) (shuffle), \(1.1\times10^{-131}\) (log), \(6.5\times10^{-99}\) (polylog) | 132/99 | gengine6.py, gtab_build.log |
| 12 | gtab assembly vs engine \(P^{(5)}/u\) | \(1.13\times10^{-96}\) | 96 | dod5d_assembly.py |
| 13 | all 30 letter differences \(S\)-unit | exact, asserted | — | qs3_symbol.py, step_55 |
| 14 | controls: Li₅ (13 pts), 300 shuffles | pair to 0, exact | — | step_72/73 |
| 15 | control: M4c \(W_5\) → (84,−38√2) | exact | — | obstruction32.py, step_79 |
| 16 | control: octagon → (−1485+1035√2, 900−540√2) | exact (and local = 0, single = 2× single-p) | — | step_79 |
| 17 | G-symbol identities (Li₂, Li₅, shuffle) | exact | — | step_73 |
| 18 | total symbol, two independent assemblies | 243/243 entries identical | — | step_90/97 |
| 19 | **obstruction table nonzero (24/24); local = 0; mirror equalities; block sum = total** | exact | — | obstruction32.py, step_81–92 |
| 20 | Li_{3,2}/Li_{4,1} correction, 24 terms | residual exactly 0 (48 eqs, rank 48) | — | step_102–109 |

## Files

- `verify_dod5d.py` — deliverable verifier (< 5 min), see §below.
- `dod_engine.py` (value engine), `dod5d_exact.py` (exact ℚ(√3) expansion), `dod5d_grewrite.py`
  (exact G-rewrite), `gengine6.py` (six-letter G-table), `qs3_symbol.py` (constant-retaining
  symbol engine over ℚ(√3) S-units), `obstruction32.py` (general antisymmetrised (3,2) cobracket
  pairings for dim V ≥ 2, any exact coefficient ring), `dod5d_assembly.py`, `itab_build2.py`,
  `dod5d_symbol_assembly.py`, `gsym_build.py`, `gtab_build.py`.
- Data: `p5_run120.txt`, `p5_run170.txt`, `p5_170.pkl`, `E5_terms.pkl`, `E5_dod_expanded.pkl`,
  `itab140.pkl`, `rewrite_exact.pkl`, `gtab6_d130_N320.pkl`, `gsym_tab.pkl`, `dod5d_symbols.pkl`,
  `mpl_cobrackets.pkl`, `li32_solution_lines.pkl`.

**Labels:** PROVED — exact expansion/rewrite/factorisations/symbols/obstruction tables/controls/
Li₃,₂ correction (symbol level; the number-level statement needs the period conjecture).
NUMERICAL — all value-level checks (engine, quadratures, G-table, assemblies), with residuals as printed.
'''
with open("result.md","w") as f: f.write(result_md)
print("result.md written:", len(result_md), "chars")
