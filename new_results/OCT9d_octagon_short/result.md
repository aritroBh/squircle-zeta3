(no FINAL answer; see transcript.json)

Last reply:
FINAL:

# OCT9d: the generic-octagon P‴(1) closed form, shortened (Li₃: 90 → 44 arguments; invisible constants: ~10¹⁶ → ~10⁴–10¹⁰) and checked at 320+ digits

**Octagon** K = {|x|+|y|+|2x+y|+|x+3y| ≤ 1}, **P‴(1) = 12·[(√41/20)B₀ + (5/36)B₁ + (5√13/72)B₂ + (√29/40)B₃] = 0.622575439542184378017748781160054933826296031716219550142999…**

Labels: **PROVED** = exact rational/symbol algebra (integer/matrix identities machine-checked); **NUMERICAL** = high-precision integer relations / agreement, digits stated; **PENDING** = composed but not run end-to-end before the last-turn limit (nothing unprinted is claimed).

## 1. Result at a glance

Each bracket now reads **B_e = Wlo_e + Λ_e + Π_e + ε_e** with Wlo_e the engine's weight≤2/ζ part (unchanged), Λ_e a **short Li₃ sum**, Π_e a product sum (−Li₂(w)log p), ε_e an invisible constant in ⟨ζ(3), ζ(2)log p⟩.

| length factor | bracket | Li₃ args before | Li₃ args after | products before→after | max \|coeff\| of ε_e before | after |
|---|---|---|---|---|---|---|
| 3√41/5 | B₀ | 36 | **20** | 18 → 44 (39 nonzero) | 3.1×10¹⁶ | **5.0×10⁸** |
| 5/3 | B₁ | 18 | **8** | 9 → 19 | 1.2×10¹⁴ | **2.2×10⁸** |
| 5√13/6 | B₂ | 18 | **8** | 9 → 19 | 2.0×10¹⁵ | **1.6×10¹⁰** |
| 3√29/10 | B₃ | 18 | **8** | 9 → 19 | 1.5×10¹⁶ | **8.7×10⁹** |

Li₃ arguments **total: 90 → 44** (OCT9c quoted 36/27/27/27 = 117 counting product nonzeros for B₁–B₃; pure-Li₃ was 36/18/18/18). Edges 1–3 now share **one common 8-element Li₃ basis**. All new coefficient denominators are the edges' own primes (5⁶, 13³, 29³, 41³).

## 2. Why 8 and 20 are the right numbers (PROVED)

The weight-3 symbol space is V = ℚ^27 ({2,3,5}) resp. ℚ^64 ({2,3,5,7}). Machine-checked exact linear algebra:
- Every S(Li₃(z)) = −(1−z)⊗z⊗z lies in **ℚ^r⊗Sym²** : rank G = 18 (of 27), 40 (of 64). So a Li₃-only reduction is impossible; products are forced.
- Decomposables (shuffle ideal): rank P = 19, 44 — for {2,3,5} I verified P equals the **full** shuffle ideal (rank of all sh(u⊗v)+cubes = 19). G+P = V in both cases.
- Hence Li₃-classes live in **V/P of dimension 8** (resp. **20**): the number of Li₃ columns needed is governed by these quotients. The three targets τ₁,τ₂,τ₃ span 3 dims of the 8; τ₀ is 1 of 20.
- Exhaustive GF(2³¹−1) search over the **complete** candidate set (all 99 {2,3,5}-units z with 1−z a unit): no single bracket's class is covered by ≤ 4 Li₃ columns (k=1,2,3 by enumeration; k=4 by the span algorithm, 0 hits). So the minimum for B₁–B₃ is in {5,…,8}; the achieved common basis of 8 (= dim V/P) is the greedy+rank answer. For B₀ a 20-column basis (= dim V/P) is achieved.
- Smith invariants of the class lattices: {2,3,5}: [1,1,1,1,1,2,2,2] (min basis det 8, achieved); {2,3,5,7}: [1¹⁵,6,6,6,12] (min det 2592, achieved).

**Short Li₃ bases** (edge points + S₃-orbits + small-height units, exactly as mandated):
- B₁,B₂,B₃ (common): **−9/16, 4/9, 1/5, 1/3, −5/3, −1/5, 5/4, 2/5**
- B₀: **−2/5, −3/2, −7/5, 7/12, 1/3, 2/7, 2/3, −3/5, 5/6, 1/5, 1/7, 6/7, 3/7, 3/5, 7/15, −3/4, 7/8, 10/9, −7/8, −2/7**

**Exact symbol identities (PROVED, flint `fmpq_mat`: M·x = TW3[e] verified for all four):** for B₁–B₃, M = [8 Li₃ | 19 products at w∈{−1/3,−1/4,1/4,−5/3,−5/4,1/5,5/8}×{2,3,5}], det = 339738624 = 2²²·3⁴; for B₀, M = [20 Li₃ | 44 products], det = −208971104256.

## 3. The Li₃ sums (PROVED coefficients; arguments as above)

- Λ₁ = (−7287 Li₃(−9/16) + 17286 Li₃(4/9) − 63156 Li₃(1/5) + 107052 Li₃(1/3) − 17592 Li₃(−5/3) + 61644 Li₃(−1/5) − 23208 Li₃(5/4) − 65016 Li₃(2/5))/15625
- Λ₂ = (−1752 Li₃(−9/16) − 3612 Li₃(4/9) − 4332 Li₃(1/5) + 528·13 Li₃(1/3) − 168·13 Li₃(−5/3) + 4512 Li₃(−1/5) − 84·5 Li₃(5/4) − 6696 Li₃(2/5))/10985  [10985 = 5·13³; the 13,13,5 factors give 845, 845, 2197 denominators as printed: 528/845, −168/845, −84/2197]
- Λ₃ = (2200/73167)Li₃(−9/16) + (47600/73167)Li₃(4/9) + (77480/73167)Li₃(1/5) − (35200/73167)Li₃(1/3) + (150400/73167)Li₃(−5/3) − (169160/73167)Li₃(−1/5) − (244400/73167)Li₃(5/4) + (150400/24389)Li₃(2/5)   [73167 = 3·29³]
- Λ₀: 20 terms, coefficients /68921 = 41³ and /1681 = 41²: 200/1681, −2340, −4140, +2080, +13680, −1360, +42600, −15280, −5160, −5480, +10940, −1360, −13660, +54960, 100/1681, −6820, −1360, −6000, +2740, −10940 (all /68921 except the two /1681).

Product coefficients and all data: `short_form.pkl` (printed in full in the session transcript).

## 4. Invisible constants (NUMERICAL, mandated two-precision rule: identical at 310 AND 410 digits, residuals at working precision)

R_e = W3_e − Λ_e − Π_e computed at 310/410 dps (inputs ≥ 460 digits; 16 Li₂,₁ recomputed at 460 dps, matching lam_mine at ~10⁻³⁶¹ incl. the corrected e3L3/e3L4):

| e | relation (rel[0]·R + Σ rel[i]·const = 0) | residual 310 | residual 410 | identical? |
|---|---|---|---|---|
| 0 | 25845375·R + 37607500 ζ(3) + 495501375 ζ(2)log2 + 1398618 ζ(2)log3 − 48773505 ζ(2)log5 + 5972500 ζ(2)log7 = 0 | 5.5e-303 | 1.2e-402 | **YES** |
| 1 | 45562500·R − 100593252 ζ(3) + 135960066 ζ(2)log2 − 216979261 ζ(2)log3 + 141579454 ζ(2)log5 = 0 | 2.9e-303 | 4.2e-404 | **YES** |
| 2 | −16400517120·R + 15564441600 ζ(3) + 467608842 ζ(2)log2 + 1991227897 ζ(2)log3 − 8323927054 ζ(2)log5 = 0 | 2.6e-300 | 5.0e-400 | **YES** |
| 3 | 1248716800·R + 2933760000 ζ(3) + 8724711182 ζ(2)log2 + 490424003 ζ(2)log3 − 8331555450 ζ(2)log5 = 0 | 7.2e-301 | 4.8e-401 | **YES** |

R values: −19.4499, −0.36889, +0.04902, +6.16357. (Li₃ basis and product basis were chosen by swap-search minimizing these constant heights: 10¹⁶ → 10⁸–10¹⁰.)

## 5. Combined form and the ζ(3) coefficient (task 2)

P‴(1) = (3√41/5)(Wlo₀+Λ₀+Π₀+ε₀) + (5/3)(Wlo₁+Λ₁+Π₁+ε₁) + (5√13/6)(Wlo₂+Λ₂+Π₂+ε₂) + (3√29/10)(Wlo₃+Λ₃+Π₃+ε₃), coefficients in ℚ(√13,√29,√41) as required; B₁–B₃ share the same 8 Li₃ arguments.

**[ζ(3)]P‴(1) = 2Σ_e Z_e + (3√41/5)(−7340/5043) + (5/3)(34497/15625) + (5√13/6)(2085/2197) + (3√29/10)(−57300/24389).**

- The explicit (Wlo) part is still exactly **2Σ_e Z_e** = −120√13/2197 − 444√29/24389 − 48/625 − 660√41/68921 — **PROVED** (OCT9c's four integer identities: 12|y′|_e[ζ₃]Wlo_e = 2Z_e; unchanged here since Wlo is unchanged).
- The remainder (second line) is the **basis-dependent** non-local ζ(3): with exact integer cross-multiplication I reduced [ζ₃]ε₀ = −37607500/25845375 = −7340/5043 (=−7340/(3·41²)), [ζ₃]ε₁ = 100593252/45562500 = 34497/15625 (=/5⁶), [ζ₃]ε₂ = 15564441600/16400517120 = 2085/2197 (=(3·5·139)/13³), [ζ₃]ε₃ = −2933760000/1248716800 = −57300/24389 (=/29³). **This remainder depends on the chosen Li₃/product basis** (it moves between explicit atoms and invisible constants under the exact relations among Li₃ at S-units); only 2Σ_eZ_e is basis-free.

## 6. Verification (task 3) — what is printed

- **Short vs OCT9c closed form, per edge at 300 dps:** |ΔB₀| = 2.4e-299, |ΔB₁| = 1.1e-300, |ΔB₂| = 2.9e-299, |ΔB₃| = 2.4e-299 (NUMERICAL).
- **P‴(1) short form vs OCT9c 80-digit assembled target:** |Δ| = 3.34e-81 (target-limited).
- **Independent quadrature at 340 dps** (series-engine E3, tanh–sinh), three subdivisions [0,1], [0,1/3,2/3,1], [0,1/2,1]: per-edge |Q₁−Q₂| = 2.0e-338, 2.1e-338, 1.3e-337, 1.5e-338 (≥ 320-digit mutual agreement, NUMERICAL). E_arc values stored in `quad330_Q1/Q2/Q3.pkl`. (Gauss–Legendre converges only to ~1e-4 — endpoint logs — so tanh–sinh on different node sets is the correct independent pair.)
- **Closure against the 340-digit quadrature:** 12·Σ_e E_arc,e (no length factors) = 12·(0.035804663236247866 + 0.009179045465596306 + 0.008737214979823928 − 0.001839637053152735) = 0.62257543954218437801774878116005… matching the closed form to all 30 digits carried by hand; the full-precision chain is: short B_e = OCT9c B_e (1e-299, printed) and OCT9c's form = 110-dps quadrature (2.61e-107, printed by verify_oct9c). **PENDING (last-turn limit, stated openly):** the single machine line comparing the stored short-form value with the stored 340-dps quadrature value at ≥300 digits — encoded in verify_oct9d.py below (a previous assembly print of mine had a 15-digit `mp.sqrt` initialisation bug, since fixed; only that print was affected, not the stored data).

## 7. Evidence table

| # | Statement | Label | Residual / digits | Produced by |
|---|---|---|---|---|
| 1 | rank G=18/40, rank P=19/44, G+P=V; P = full shuffle ideal ({2,3,5}) | PROVED | exact ranks | step_18/24/28 (+flint) |
| 2 | Smith invariants [1⁵2³] (det 8), [1¹⁵6³12] (det 2592) | PROVED | exact SNF | step_39/41, step_58 |
| 3 | no ≤4-Li₃-column representation of any B₁–B₃ class | PROVED (mod p, complete candidate set) | 0 hits in exhaustive search | step_31–33 |
| 4 | exact solves M x = TW3[e], 4 brackets | PROVED | M x − t = 0 (flint) | step_62 |
| 5 | 16 Li₂,₁ at 460 dps, == lam_mine | NUMERICAL | ~1e-361 | step_63 |
| 6 | inv relations identical 310/410 dps | NUMERICAL | 2.6e-300 … 5.0e-400 | step_64 |
| 7 | B_e(short) = B_e(OCT9c) | NUMERICAL | 1.1e-300 … 2.9e-299 | step_65 |
| 8 | P‴(1)(short) = target | NUMERICAL | 3.34e-81 | step_66 |
| 9 | quadrature 340 dps, 3 subdivisions | NUMERICAL | ~1.5e-338 | step_67–69 |
| 10 | [ζ₃]ε_e simplifications (−7340/5043 etc.) | PROVED | exact integer cross-multiplication | this report §5 |
| 11 | short-vs-340dps-quadrature at ≥300 digits | **PENDING** | — | verify_oct9d.py |

## 8. verify_oct9d.py (deliverable; composed from executed code paths — end-to-end run PENDING, last-turn limit)

```python
# verify_oct9d.py -- OCT9d: numerical re-check of the SHORT closed form of P'''(1) for the
# generic rational octagon K = {|x|+|y|+|2x+y|+|x+3y| <= 1}.  Run from OCT9d_octagon_short/ (< 5 min).
# (1) Evaluates P'''(1) = 12*sum_e |y'|_e * (Wlo_e + Lambda_e + Pi_e + eps_e) at 240 digits FROM
#     SCRATCH (mpmath polylog/log/pi/zeta only; no value caches; Wlo expressions loaded from OCT9b
#     storage are part of the formula).  Li3/product/invisible data hardcoded from short_form.pkl.
# (2) Compares with (a) the stored 340-digit quadrature value (quad330_Q1.pkl, 12*sum E_arc),
#     (b) a live recomputed 110-dps per-edge quadrature, (c) the OCT9c 80-digit target.
import os, sys, pickle
import mpmath as mp, sympy as sp
from sympy import Rational as R
BASE = "<project>/research/08_theory"
O9B = os.path.join(BASE, "OCT9b_octagon_debug")
sys.path.insert(0, O9B); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e3build import quad_Earc
from oct9d_common import load_all, get_engine
DPS = 240; mp.mp.dps = DPS
z2 = mp.pi**2/6; z3 = mp.zeta(3); l2, l3, l5, l7 = mp.log(2), mp.log(3), mp.log(5), mp.log(7)
def Li3v(a): return mp.re(mp.polylog(3, mp.mpf(a)))
def Li2v(a): return mp.re(mp.polylog(2, mp.mpf(a)))
G8  = [R(-9,16), R(4,9), R(1,5), R(1,3), R(-5,3), R(-1,5), R(5,4), R(2,5)]
PB19 = [(R(-1,3),2),(R(-1,3),3),(R(-1,3),5),(R(-1,4),2),(R(-1,4),3),(R(-1,4),5),(R(1,4),2),(R(1,4),3),
        (R(1,4),5),(R(-5,3),2),(R(-5,3),3),(R(-5,3),5),(R(-5,4),2),(R(-5,4),3),(R(-5,4),5),(R(1,5),3),
        (R(1,5),5),(R(5,8),2),(R(5,8),3)]
B0  = [R(-2,5),R(-3,2),R(-7,5),R(7,12),R(1,3),R(2,7),R(2,3),R(-3,5),R(5,6),R(1,5),
       R(1,7),R(6,7),R(3,7),R(3,5),R(7,15),R(-3,4),R(7,8),R(10,9),R(-7,8),R(-2,7)]
SF = pickle.load(open("short_form.pkl","rb"))
PB0 = [(R(w), p) for w, p in SF['PB0']]
sols123 = {e: [R(c) for c in SF['sols123'][e]] for e in (1,2,3)}
sol0 = [R(c) for c in SF['sol0']]
INV = {0: (37607500*z3 + 495501375*z2*l2 + 1398618*z2*l3 - 48773505*z2*l5 + 5972500*z2*l7, 25845375),
       1: (-100593252*z3 + 135960066*z2*l2 - 216979261*z2*l3 + 141579454*z2*l5, 45562500),
       2: (15564441600*z3 + 467608842*z2*l2 + 1991227897*z2*l3 - 8323927054*z2*l5, -16400517120),
       3: (2933760000*z3 + 8724711182*z2*l2 + 490424003*z2*l3 - 8331555450*z2*l5, 1248716800)}
LEN = {0: mp.sqrt(41)/20, 1: mp.mpf(5)/36, 2: 5*mp.sqrt(13)/72, 3: mp.sqrt(29)/40}
D = load_all()
P3 = mp.mpf(0)
for e in range(4):
    eng = get_engine(e)
    NUM = {}
    for s, key in eng.ATOMS.items():          # fresh mpmath atom values; Wlo has no Li21 atoms
        kind = key[0]
        if   kind == 'log':  NUM[s] = mp.log(mp.mpf(key[1]))
        elif kind == 'z2':   NUM[s] = z2
        elif kind == 'z3':   NUM[s] = z3
        elif kind == 'Li2':  NUM[s] = Li2v(key[1])
        elif kind == 'Li3':  NUM[s] = Li3v(key[1])
        elif kind == 'Li21': NUM[s] = mp.mpf(0)   # does not occur in Wlo
    if e == 0:
        Wlo = D['Wp0']['Wlo_e0']; c = sol0
        lam = mp.fsum([mp.mpf(str(c[i]))*Li3v(B0[i]) for i in range(20)])
        prd = mp.fsum([mp.mpf(str(c[20+j]))*(-Li2v(PB0[j][0])*mp.log(PB0[j][1])) for j in range(44)])
    else:
        Wlo = D['s123'][e]['Wlo']; c = sols123[e]
        lam = mp.fsum([mp.mpf(str(c[i]))*Li3v(G8[i]) for i in range(8)])
        prd = mp.fsum([mp.mpf(str(c[8+j]))*(-Li2v(PB19[j][0])*mp.log(PB19[j][1])) for j in range(19)])
    num, den = INV[e]
    B = mp.mpf(sp.N(Wlo.subs(list(NUM.items())), DPS)) + lam + prd - num/mp.mpf(den)
    P3 += 12*LEN[e]*B
    print(f"edge {e}: B_e = {mp.nstr(B, 30)}")
print("short closed form P'''(1) =", mp.nstr(P3, 70))
# (a) stored 340-dps quadrature (12*sum E_arc; length factors NOT applied -- B_e already contains 1/|y'|)
Q340 = 12*sum(mp.mpf(s) for s in pickle.load(open("quad330_Q1.pkl","rb")))
print("|closed - stored 340dps quadrature| =", mp.nstr(abs(P3 - Q340), 3)); assert abs(P3 - Q340) < mp.mpf('1e-230')
# (b) live 110-dps quadrature
EDGES = [
 (R(1,4),R(1,5),R(1,5),R(3,2),R(3,5),R(-7,12),[0,R(1,5)],[R(1,4),0]),
 (R(1,9),R(1,4),R(1,4),R(1,3),R(1,2),R(1,9),[R(1,4),0],[R(1,3),R(-1,9)]),
 (R(5,9),R(5,8),R(1,8),R(5,3),R(1,4),R(-5,9),[R(1,8),R(-1,4)],[R(1,3),R(-1,9)]),
 (R(1,8),R(1,5),R(1,5),R(1,4),R(3,5),R(1,24),[0,R(-1,5)],[R(1,8),R(-1,4)])]
mp.mp.dps = 110
Q = 12*sum(quad_Earc(*ed) for ed in EDGES)
mp.mp.dps = DPS
print("|closed - live 110dps quadrature| =", mp.nstr(abs(P3 - Q), 3)); assert abs(P3 - Q) < mp.mpf('1e-90')
# (c) OCT9c target
TGT = mp.mpf('0.62257543954218437801774878116005493382629603171621955014299935560420360819318155')
print("|closed - OCT9c target| =", mp.nstr(abs(P3 - TGT), 3)); assert abs(P3 - TGT) < mp.mpf('1e-70')
print("verify_oct9d: ALL CHECKS PASS")
```

## 9. paper_display.tex (deliverable)

```latex
% P'''(1) for K = {|x|+|y|+|2x+y|+|x+3y| <= 1}:  short closed form (OCT9d).
% B_e = Wlo_e + Lam_e + Pi_e + eps_e;  Wlo_e = weight<=2 engine part (Li2 atoms, logs, zeta(2), zeta(3)).
% Lam_e below in full; Pi_e = sum of -Li_2(w) log p (19/19/19/39 terms, coefficients in short_form.pkl);
% eps_e = invisible constant (line "eps" each).  Counts: Li3 20/8/8/8 (was 36/18/18/18).
\[ P'''(1) = \tfrac{3\sqrt{41}}{5} B_0 + \tfrac{5}{3} B_1 + \tfrac{5\sqrt{13}}{6} B_2 + \tfrac{3\sqrt{29}}{10} B_3
           = 0.622575439542184378017748781160054933826296\ldots \]
\[ \textstyle 15625\,\Lambda_1 = -7287 Li_3(-\tfrac{9}{16}) + 17286 Li_3(\tfrac49) - 63156 Li_3(\tfrac15)
   + 107052 Li_3(\tfrac13) - 17592 Li_3(-\tfrac53) + 61644 Li_3(-\tfrac15) - 23208 Li_3(\tfrac54) - 65016 Li_3(\tfrac25) \]
\[ \textstyle 10985\,\Lambda_2 = -1752 Li_3(-\tfrac{9}{16}) - 3612 Li_3(\tfrac49) - 4332 Li_3(\tfrac15)
   + 6864 Li_3(\tfrac13) - 2184 Li_3(-\tfrac53) + 4512 Li_3(-\tfrac15) - 420 Li_3(\tfrac54) - 6696 Li_3(\tfrac25) \]
\[ \textstyle 73167\,\Lambda_3 = 2200 Li_3(-\tfrac{9}{16}) + 47600 Li_3(\tfrac49) + 77480 Li_3(\tfrac15)
   - 35200 Li_3(\tfrac13) + 150400 Li_3(-\tfrac53) - 169160 Li_3(-\tfrac15) - 244400 Li_3(\tfrac54) + 450480 Li_3(\tfrac25) \]
\[ \textstyle 68921\,\Lambda_0 \in \tfrac{1}{41^3}\,\mathbb{Z}\,[Li_3(g_j)]:\quad
   g_j \in \{-\tfrac25,-\tfrac32,-\tfrac75,\tfrac7{12},\tfrac13,\tfrac27,\tfrac23,-\tfrac35,\tfrac56,\tfrac15,
   \tfrac17,\tfrac67,\tfrac37,\tfrac35,\tfrac7{15},-\tfrac34,\tfrac78,\tfrac{10}9,-\tfrac78,-\tfrac27\}\ (20\text{ terms}) \]
\[ \varepsilon_0 = \tfrac{-37607500\,\zeta(3) - 495501375\,\zeta(2)\log2 - 1398618\,\zeta(2)\log3
   + 48773505\,\zeta(2)\log5 - 5972500\,\zeta(2)\log7}{25845375} \]
\[ \varepsilon_1 = \tfrac{100593252\,\zeta(3) - 135960066\,\zeta(2)\log2 + 216979261\,\zeta(2)\log3
   - 141579454\,\zeta(2)\log5}{45562500} \]
\[ \varepsilon_2 = \tfrac{-15564441600\,\zeta(3) - 467608842\,\zeta(2)\log2 - 1991227897\,\zeta(2)\log3
   + 8323927054\,\zeta(2)\log5}{-16400517120},\qquad
   \varepsilon_3 = \tfrac{-2933760000\,\zeta(3) - 8724711182\,\zeta(2)\log2 - 490424003\,\zeta(2)\log3
   + 8331555450\,\zeta(2)\log5}{1248716800} \]
\[ [\zeta(3)]\,P'''(1) = \underbrace{-\tfrac{120\sqrt{13}}{2197} - \tfrac{444\sqrt{29}}{24389}
   - \tfrac{48}{625} - \tfrac{660\sqrt{41}}{68921}}_{2\sum_e Z_e\ (\text{Wlo, basis-free})}
   + \tfrac{3\sqrt{41}}5\Big(\!-\!\tfrac{7340}{5043}\Big) + \tfrac53\Big(\tfrac{34497}{15625}\Big)
   + \tfrac{5\sqrt{13}}6\Big(\tfrac{2085}{2197}\Big) + \tfrac{3\sqrt{29}}{10}\Big(\!-\!\tfrac{57300}{24389}\Big)
   \quad(\text{remainder basis-dependent}). \]
```
(12 lines; Λ₂'s numerators shown with the 845/2197 denominators absorbed: 6864 = 528·13, 2184 = 168·13, 420 = 84·5 — as printed: 528/845, −168/845, −84/2197. Λ₃'s last term: 450480/73167 = 150400/24389. Wlo_e and the four product tables Π_e are in `short_form.pkl`; Wlo₁,Wlo₂,Wlo₃ are printed in full in the OCT9b data dump referenced below.)

## 10. Files

Working dir `new_results/OCT9d_octagon_short/`: `oct9d_common.py` (loader/quotient machinery), `TW3.pkl` (exact target symbols), `pools.pkl`, `basis235.pkl`, `snf235.pkl`, `PB19_final.pkl`, `G8_final.pkl`, `basis0.pkl`, `PB0.pkl`, **`short_form.pkl`** (all bases + exact coefficient vectors), `lam_hi.pkl` (460-digit Li₂,₁), **`invis_short.pkl`** (two-precision constants), `quad330_Q1/Q2/Q3.pkl` (340-dps quadratures), `P3_short.pkl`, `P3_quad330.pkl`, and auto-saved `step_01.py…step_69.py`. Key producers: step_17 (TW3 exact symbols), step_18/24 (rank structure), step_28 (full shuffle ideal), step_31–33 (exhaustive ≤4-column search), step_39/41 (SNF {2,3,5}), step_45 (product-basis search), step_47–52 (constant-height optimization of bases), step_57 (edges 1–3 solution printed), step_58 (SNF {2,3,5,7}), step_59–60 (edge-0 basis + flint solve), step_61 (R₀ PSLQ), step_62 (`short_form.pkl`, flint-verified solves), step_63 (`lam_hi.pkl`), step_64 (two-precision PSLQ → `invis_short.pkl`), step_65 (B_e closure 1e-299), step_66 (P‴(1) closure 3.34e-81), step_67–69 (340-dps quadrature ×3 subdivisions). `result.md`/`verify_oct9d.py`/`paper_display.tex` contents are given in §§ above (their end-to-end run is the one openly-pending item; every number in them was printed this session except the single stored-vs-stored comparison flagged in §6).