# RUL1 — exact level-2 reductions at weight ≤ 5: three "numerical" ingredients turned into proofs

**Deliverable directory:** `new_results/RUL1_level2_exact_reductions` · **`verify_rul1.py`: 41/41 checks PASS in 235 s (< 5 min)** — rebuilds the relations, reduces every target with certificates, checks each numerically to ≥100 digits.

## 1. The exact reduction engine (all relations proved, never PSLQ)

Alternating MZVs (level 2) as iterated integrals $G(w;1)$, words over $\{0,1,-1\}$ (letter $a$ = form $dt/(t-a)$), with $\zeta(S)=(-1)^dG(w(S);1)$, $e_i=\varepsilon_1\cdots\varepsilon_i$. Relation families, each an instance of a **proved** theorem applied with rational arithmetic:

- **DS** — finite double shuffle (shuffle of convergent words minus transported stuffle; from absolute convergence);
- **DUAL** — the involution $t\mapsto(1-t)/(1+t)$ on convergent words ($0\mapsto(1)-(-1)$, $1\mapsto(0)-(-1)$, $-1\mapsto-(-1)$, $G(w)=(-1)^{|w|}G(\tau^*w^{rev})$; a change of variables). No regularization/IKZ theorem is needed.

Exact ranks (Fraction Gauss elimination with history): codimensions of convergent-word spaces = **1, 2, 3, 5, 8** at weights 1–5 = the Deligne–Goncharov dimensions:

| weight | conv. words | rels (DS+DUAL) | rank | codim |
|---|---|---|---|---|
| 1 | 1 | 1 | 0 | 1 |
| 2 | 4 | 5 | 2 | 2 |
| 3 | 12 | 16 | 9 | 3 |
| 4 | 36 | 58 | 31 | 5 |
| 5 | 108 | 192 | 100 | 8 |

Named basis (the task's): $\log2$; $\pi^2,\log^22$; $\zeta(3),\pi^2\log2,\log^32$; $\mathrm{Li}_4(\tfrac12),\zeta(3)\log2,\pi^4,\pi^2\log^22,\log^42$; $\mathrm{Li}_5(\tfrac12),\mathrm{Li}_4(\tfrac12)\log2,\zeta(5),\zeta(3)\pi^2,\zeta(3)\log^22,\pi^4\log2,\pi^2\log^32,\log^52$, with $\mathrm{Li}_n(1/2)$ as explicit $2^{n-1}$-word combos via $t=s/(1+s)$ (70-digit checked). Named-basis determinants $12, 24, -11520/7, 1658880/91 \neq 0$ ⇒ reductions unique. **Every target identity carries a certificate** $T-\sum c_i\mathrm{named}_i=\sum_j h_j R_j$ over proved relation instances, verified by exact rational arithmetic.

## 2. Controls (PASS)

$\zeta(\bar2,1)=\frac18\zeta(3)$ (certificate) ⇒ $\sum(-1)^{n+1}H_n/n^2=\eta(3)-\zeta(\bar2,1)=\frac58\zeta(3)$; $\mathrm{Li}_4(1/2)$ primitive; $\zeta(\bar3,\bar1)=-2\mathrm{Li}_4(\frac12)+\frac{\pi^4}{180}+\frac{\pi^2\log^22}{12}-\frac{\log^42}{12}$ (certificate; 121-digit agreement); dimension counts.

## 3. Target 1 — the ruler $\pi_p$ (thm:pip): $N_4,N_5$ **PROVED**

The "five classical alternating Euler sums" are exactly $J_{a,k}=a!\,k!\,G(0^{k+1},(-1)^a;1)=(-1)^aa!k!\,\zeta(\overline{k+2},1^{a-1})$:

| J | closed form (= D12's PSLQ value) | cert | numeric residual |
|---|---|---|---|
| $J_{3,0}$ | $-\frac{\log^42}4+\frac{\pi^2\log^22}4+\frac{\pi^4}{15}-\frac{21}4\zeta(3)\log2-6\mathrm{Li}_4$ | ✓ | 7.5e-121 |
| $J_{2,1}$ | $\frac{\log^42}6-\frac{\pi^2\log^22}6-\frac{\pi^4}{24}+\frac72\zeta(3)\log2+4\mathrm{Li}_4$ | ✓ | 4.8e-121 |
| $J_{4,0}$ | $-\frac45\log^52+\frac23\pi^2\log^32-\frac{21}2\zeta(3)\log^22+24\zeta(5)-24\mathrm{Li}_4\log2-24\mathrm{Li}_5$ | ✓ | 3.9e-120 |
| $J_{3,1}$ | $\frac25\log^52-\frac13\pi^2\log^32+\frac{21}4\zeta(3)\log^22-\frac12\zeta(3)\pi^2-\frac{99}{16}\zeta(5)+12\mathrm{Li}_4\log2+12\mathrm{Li}_5$ | ✓ | 2.0e-120 |
| $J_{2,2}$ | $\frac13\zeta(3)\pi^2-\frac{29}8\zeta(5)$ | ✓ | 2.2e-121 |

$\frac1{12}J_{3,0}+\frac18J_{2,1}=\frac{\pi^4}{2880}$ and $\frac1{48}(J_{4,0}+2J_{3,1})=-\frac{\zeta(3)\pi^2}{48}+\frac{31\zeta(5)}{128}$ hold **exactly**. Chain to $N_4,N_5$ with **no numerical step**: D12's exact formula $N_k=\frac{(-\log2)^{k+1}}{(k+1)!}+\sum_{m=0}^{k-1}\int[h^m]\Psi$ rebuilt rationally; only $J_{1,k}=k!\,\eta(k+2)$ (elementary), $J_{2,0}=\zeta(3)/4$ (control) and the five proved J's enter ($\mathrm{Li}_4,\mathrm{Li}_5$ cancel identically), yielding exactly the D12 closed forms $N_4,N_5$; independent quadrature assembly agrees to 5.6e-132 / 1.7e-131.

## 4. Target 2 — non-local $\zeta(5)$ of the hexagon edge $r=2$: **PROVED**

All 13 HPLs at $x=-1$ ($H_w(-1)=(-1)^{\#1}G(w[1{\mapsto}{-1}];1)$) reduced with certificates, matching K1g's former PSLQ values, 120-digit confirmed (e.g. $H_{00111}=\frac{33}{32}\zeta(5)+\frac1{12}\zeta(3)\pi^2-\frac78\zeta(3)\log^22+\frac1{18}\pi^2\log^32-\frac1{15}\log^52-2\mathrm{Li}_5-2\mathrm{Li}_4\log2$). **Exact fold** of K1g's 307-monomial assembly (543 exact reductions):
$$[\zeta(5)]E_{arc}=-\tfrac{20187}{2000}\sqrt5,\quad [\zeta(2)\zeta(3)]E_{arc}=-\tfrac{858}{125}\sqrt5 \ \textbf{(PROVED)};\quad N^{(5)}=-\tfrac{711}{400}\sqrt5,\ N^{(2,3)}=-\tfrac{1239}{250}\sqrt5$$
(the far kink lines, after subtracting the exact local law $-\frac{2079}{250}\sqrt5$, $-\frac{477}{250}\sqrt5$). The conv-B fold ($\ell_j=\log(A_j/|m_j|)=(\log2,-\frac12\log2,0)$, log 3 never appears) equals $\sqrt5\,X$ with E5e's full 20-constant closed form **exactly**, and matches the independent 337-digit arc value to 2.96e-123.

## 5. Target 3 — regular hexagon ($r=1$): **PROVED**

Paper's hexagon ($w_j=1$, edge $1/\sqrt3$): edge data $A_0=A_1=B=\frac12,\kappa=0,L_2=\frac13,q_0=-\frac16,V_{02}=\frac13$, logs $-\log2$. At $r=1$ only local terms survive; the assembly gives $P^{(5)}(1)=6E_{arc}$ equal **exactly** to the paper's 15-term closed form (128/√3·P⁽⁵⁾(1) = 30720 log2 − 4920π² − … + 256 log⁵2), with $[\zeta(5)]=\frac{2835\sqrt3}{32}$ and $[\zeta(2)\zeta(3)]=\frac{1395\sqrt3}{32}$ exactly; matches D17's independent 100-digit value to 4.98e-98.

## 6. Evidence table (printed by the listed code)

| Check | Result | Residual | Digits | Code |
|---|---|---|---|---|
| relation codims = DG | 1,2,3,5,8 | exact | exact | verify_rul1, rul1_step02 |
| controls (ζ(\bar2,1), ζ(\bar3,bar1), Li_4(1/2), J_{2,0}) | §2 | certs exact; 1.3e-121 | 120 | rul1_step03 |
| five J_{a,k} | §3 | certs exact; ≤9.7e-120 | 120 | rul1_step04 |
| N₄, N₅ closed forms | §3 | exact vs D12 formula; 5.6e-132, 1.7e-131 | 130 | rul1_step04b |
| 13 atoms @ x=−1 | §4 | certs exact; ≤1.8e-120 | 120 | rul1_step05 |
| fold ζ5, ζ2ζ3 = −20187√5/2000, −858√5/125 | §4 | exact | exact | rul1_step05b |
| far-kink −711√5/400, −1239√5/250 | §4 | exact | exact | rul1_step05b |
| conv-B fold = √5·X_E5e (20 consts) | §4 | exact; 2.96e-123 vs 337-digit data | 120+ | rul1_step05c |
| regular hexagon 15-term form | §5 | exact; 4.98e-98 vs D17 | 98 | rul1_step06 |

## 7. Files that produced these

`rul1.py` (engine: DS+DUAL relations), `rul1la.py` (exact Gauss + history), `rul1_named.py` (basis combos), `rul1_reduce.py` (reductions + certificates), `rul1num.py` (E1 Landen / E2 quad / E3 evaluators), `rul1_step01.py` (sanity), `rul1_step02.py` (relations→`rul1_data.pkl`), `rul1_step03.py` (controls), `rul1_step04.py` (J's→`rul1_certs_J.pkl`), `rul1_step04b.py` (N₄,N₅), `rul1_step05.py` (atoms→`rul1_certs_atoms.pkl`), `rul1_step05b.py` (polar fold, ζ-coords), `rul1_step05c.py` (conv-B fold vs E5e), `rul1_step05d.py` (convention check), `rul1_step06.py` (regular hexagon), **`verify_rul1.py`** (the < 5 min verifier), `result.md`, `progress.log`.

**Labels.** PROVED: everything in §2–§5 (exact rational identities modulo proved DS/DUAL instances, certificates stored). NUMERICAL: the digit-agreement columns and the conjectural linear independence of the basis constants (which affects only the words "the" coefficients, not any identity).