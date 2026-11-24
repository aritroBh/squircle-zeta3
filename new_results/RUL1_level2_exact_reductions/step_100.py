import pathlib
md = '''# RUL1 — exact level-2 reductions at weight ≤ 5 (three "numerical" ingredients turned into proofs)

Work dir: `research/08_theory/RUL1_level2_exact_reductions`. Labels: **PROVED** (exact relations, rational arithmetic) / **NUMERICAL**.
Everything is reproduced by **`verify_rul1.py`** (rebuilds the relations, reduces every target with certificates,
checks each numerically to ≥ 100 digits): **41/41 checks PASS in 235 s**.

## 0. The engine (what "exact reduction" means here)

Alternating MZVs (level 2) as iterated integrals $G(w;1)$ over the alphabet $\\{0,1,-1\\}$
(letter $a$ = form $dt/(t-a)$), equivalently series $\\zeta(s_1,\\dots,s_d;\\varepsilon_1,\\dots,\\varepsilon_d)$,
$\\varepsilon_i=\\pm1$, with $\\zeta(S)=(-1)^d G(w(S);1)$, $e_i=\\varepsilon_1\\cdots\\varepsilon_i$ (`rul1.py`).
Relation families used, **each a proved theorem, applied with exact rational arithmetic**:
- **DS** (finite double shuffle): for convergent $u,v$, $u\\,\\unicode{x2E93}\\,v$ (word shuffle) minus the transported
  stuffle $st(u)*st(v)$ — from absolute convergence of the series/integrals.
- **DUAL** (level-2 duality): for convergent $w$, $G(w;1)=(-1)^{|w|}G(\\tau^* w^{rev};1)$ from the change of variables
  $t\\mapsto(1-t)/(1+t)$ ($0\\mapsto(1)-(-1)$, $1\\mapsto(0)-(-1)$, $-1\\mapsto-(-1)$ on letters).

No PSLQ anywhere; no regularization theorem is needed (duality supplies what RDS would, and it is elementary).
Exact rank of the relation system (sparse Fraction Gauss elimination with history, `rul1la.py`):

| weight | convergent words | relations (DS+DUAL) | rank | codimension (= Deligne–Goncharov) |
|---|---|---|---|---|
| 1 | 1 | 1 | 0 | **1** |
| 2 | 4 | 5 | 2 | **2** |
| 3 | 12 | 16 | 9 | **3** |
| 4 | 36 | 58 | 31 | **5** |
| 5 | 108 | 192 | 100 | **8** |

Named basis (`rul1_named.py`): $\\log2$, $\\pi^2$, $\\log^22$; $\\zeta(3),\\pi^2\\log2,\\log^32$;
$\\operatorname{Li}_4(\\tfrac12),\\zeta(3)\\log2,\\pi^4,\\pi^2\\log^22,\\log^42$;
$\\operatorname{Li}_5(\\tfrac12),\\operatorname{Li}_4(\\tfrac12)\\log2,\\zeta(5),\\zeta(3)\\pi^2,\\zeta(3)\\log^22,\\pi^4\\log2,\\pi^2\\log^32,\\log^52$.
$\\operatorname{Li}_n(1/2)$ enters as the explicit $2^{n-1}$-word combination via $t=s/(1+s)$ (verified to 70 digits).
Each named-basis matrix over the free columns has nonzero exact determinant
($12, 24, -11520/7, 1658880/91$ for weights 2–5), so every reduction is unique.
**Certificates:** each target $T$ comes with an exact identity $T-\\sum c_i\\,\\mathrm{named}_i=\\sum_j h_j R_j$
($R_j$ relation instances), checked by exact rational arithmetic (`rul1_reduce.py: check_certificate`).

Independent evaluators for the numeric confirmations (`rul1num.py`): E1 = Landen $x\\mapsto x/(x-1)$ to the
geometric 1/2-series (only defining series used; ≥ 120 digits); E2 = 1D quadrature of
$J_{a,k}=\\int_0^\\infty s^k\\log^a(1+e^{-s})ds$; E3 = nested-tail + parity-paired Richardson (spot checks).

## 1. Controls (all PASS)

- $\\zeta(\\bar2,1)=\\frac18\\zeta(3)$ (certificate), hence $\\sum_{n\\ge1}(-1)^{n+1}H_n/n^2=\\eta(3)-\\zeta(\\bar2,1)=\\frac34\\zeta(3)-\\frac18\\zeta(3)=\\frac58\\zeta(3)$.
- $\\operatorname{Li}_4(1/2)$ is primitive (reduces to itself, coefficient 1).
- $\\zeta(\\bar3,\\bar1)=-2\\operatorname{Li}_4(\\tfrac12)+\\frac{\\pi^4}{180}+\\frac{\\pi^2\\log^22}{12}-\\frac{\\log^42}{12}$ (certificate; matches the literature value to 121 digits).
- Dimension counts 1, 2, 3, 5, 8 (table above).

## 2. Target 1 — the ruler $\\pi_p$ (paper Thm. thm:pip): $N_4,N_5$ PROVED

The "five classical alternating Euler sums" of D12 are exactly
$J_{a,k}=a!\\,k!\\,G(0^{k+1},(-1)^a;1)$ (weight $a+k+1$): with $x=e^{-s}$,
$J_{a,k}=\\int_0^\\infty s^k\\log^a(1+e^{-s})ds=(-1)^a a!k!\\,\\zeta(\\overline{k+2},1^{a-1})$.
Exact reductions (certificates), each = the value D12 found by PSLQ:

| J | closed form | cert | |named−E1| | |named−E2| |
|---|---|---|---|---|---|
| $J_{3,0}$ | $-\\frac{\\log^42}{4}+\\frac{\\pi^2\\log^22}{4}+\\frac{\\pi^4}{15}-\\frac{21}{4}\\zeta(3)\\log2-6\\operatorname{Li}_4(\\frac12)$ | ✓ | 7.5e-121 | 7.5e-121 |
| $J_{2,1}$ | $\\frac{\\log^42}{6}-\\frac{\\pi^2\\log^22}{6}-\\frac{\\pi^4}{24}+\\frac72\\zeta(3)\\log2+4\\operatorname{Li}_4(\\frac12)$ | ✓ | 4.1e-121 | 4.8e-121 |
| $J_{4,0}$ | $-\\frac45\\log^52+\\frac23\\pi^2\\log^32-\\frac{21}{2}\\zeta(3)\\log^22+24\\zeta(5)-24\\operatorname{Li}_4(\\frac12)\\log2-24\\operatorname{Li}_5(\\frac12)$ | ✓ | 3.9e-120 | 3.9e-120 |
| $J_{3,1}$ | $\\frac25\\log^52-\\frac13\\pi^2\\log^32+\\frac{21}{4}\\zeta(3)\\log^22-\\frac12\\zeta(3)\\pi^2-\\frac{99}{16}\\zeta(5)+12\\operatorname{Li}_4(\\frac12)\\log2+12\\operatorname{Li}_5(\\frac12)$ | ✓ | 2.0e-120 | 2.0e-120 |
| $J_{2,2}$ | $\\frac13\\zeta(3)\\pi^2-\\frac{29}{8}\\zeta(5)$ | ✓ | 9.7e-122 | 2.2e-121 |

The cancellations $\\frac1{12}J_{3,0}+\\frac18J_{2,1}=\\frac{\\pi^4}{2880}$ and
$\\frac1{48}(J_{4,0}+2J_{3,1})=-\\frac{\\zeta(3)\\pi^2}{48}+\\frac{31\\zeta(5)}{128}$ hold **exactly**.
Chain to $N_4,N_5$ (**PROVED, no numerical step**): D12's exact formula
$N_k=\\frac{(-\\log2)^{k+1}}{(k+1)!}+\\sum_{m=0}^{k-1}\\int_0^\\infty[h^m]\\Psi\\,ds$ is rebuilt with rational series
arithmetic; every monomial is $s^iL^j\\log^\\ell2$ with $j\\ge1$, $i+j\\le m$, so only
$J_{1,k}=k!\\,\\eta(k+2)$ (elementary), $J_{2,0}=\\zeta(3)/4$ (control 1) and the five proved J's enter:
$$N_4=-\\frac{\\log^52}{120}+\\frac{\\pi^2}{24}\\Bigl(1-\\log2+\\frac{\\log^22}{2}\\Bigr)+\\frac{\\zeta(3)}{4}(1-\\log2)+\\frac{\\pi^4}{2880},$$
$$N_5=\\frac{\\log^62}{720}+\\frac{\\pi^2}{24}\\Bigl(1-\\log2+\\frac{\\log^22}{2}-\\frac{\\log^32}{6}\\Bigr)+\\frac{\\zeta(3)}{4}\\Bigl(1-\\log2+\\frac{\\log^22}{2}\\Bigr)+\\frac{\\pi^4}{2880}(1-\\log2)-\\frac{7\\pi^4}{1920}-\\frac{\\zeta(3)\\pi^2}{48}+\\frac{\\zeta(5)}{8}.$$
($\\operatorname{Li}_4,\\operatorname{Li}_5$ cancel identically.) Numeric: independent quadrature assembly of the
same $J$-formula agrees with the closed forms to 5.6e-132 ($N_4$) and 1.7e-131 ($N_5$).

## 3. Target 2 — non-local $\\zeta(5)$ of the hexagon edge $r=2$ ($x=-1$): PROVED

The 13 HPLs at $x=-1$: $H_w(-1)=(-1)^{\\#1}G(w[1\\mapsto-1];1)$ — alternating MZVs, reduced exactly
(certificates in `rul1_certs_atoms.pkl`), all matching K1g's former PSLQ values, confirmed to ~120 digits via E1:
$H_1=-\\log2$, $\\operatorname{Li}_2=-\\pi^2/12$, $\\operatorname{Li}_3=-3\\zeta(3)/4$, $H_{011}=\\zeta(3)/8$,
$\\operatorname{Li}_4=-7\\pi^4/720$,
$H_{0011}=2\\operatorname{Li}_4+\\frac74\\zeta(3)\\log2-\\frac{\\pi^4}{48}-\\frac{\\pi^2\\log^22}{12}+\\frac{\\log^42}{12}$,
$H_{0111}=\\operatorname{Li}_4+\\frac78\\zeta(3)\\log2-\\frac{\\pi^4}{90}-\\frac{\\pi^2\\log^22}{24}+\\frac{\\log^42}{24}$,
$\\operatorname{Li}_5=-\\frac{15}{16}\\zeta(5)$,
$H_{00011}=-\\frac{29}{32}\\zeta(5)+\\frac1{12}\\zeta(3)\\pi^2$,
$H_{00101}=\\frac{41}{32}\\zeta(5)-\\frac5{48}\\zeta(3)\\pi^2$,
$H_{00111}=\\frac{33}{32}\\zeta(5)+\\frac1{12}\\zeta(3)\\pi^2-\\frac78\\zeta(3)\\log^22+\\frac1{18}\\pi^2\\log^32-\\frac1{15}\\log^52-2\\operatorname{Li}_5-2\\operatorname{Li}_4\\log2$,
$H_{01011}=-\\frac{73}{64}\\zeta(5)-\\frac14\\zeta(3)\\pi^2+\\frac74\\zeta(3)\\log^22-\\frac19\\pi^2\\log^32+\\frac2{15}\\log^52+4\\operatorname{Li}_5+4\\operatorname{Li}_4\\log2$,
$H_{01111}=\\zeta(5)-\\frac7{16}\\zeta(3)\\log^22+\\frac1{36}\\pi^2\\log^32-\\frac1{30}\\log^52-\\operatorname{Li}_5-\\operatorname{Li}_4\\log2$.

**Exact fold** of K1g's 307-monomial assembly (`coef_acc_raw.pkl`, edge 1: $A_0=1/2,A_1=3/4,B=1/4,\\kappa=1,L_2=5,q_0=1,V_{02}=1$),
every monomial evaluated by exact shuffle products + exact reductions (543 reductions):
$$[\\zeta(5)]E_{arc}=-\\frac{20187}{2000}\\sqrt5,\\qquad [\\zeta(2)\\zeta(3)]E_{arc}=-\\frac{858}{125}\\sqrt5\\quad\\textbf{(PROVED)};$$
subtracting the local law $Z_e^{(5)}=-\\frac{2079}{250}\\sqrt5$, $Z_e^{(2,3)}=-\\frac{477}{250}\\sqrt5$ (E5_edge_law, exact):
the far kink lines add
$$N^{(5)}=-\\frac{711}{400}\\sqrt5,\\qquad N^{(2,3)}=-\\frac{1239}{250}\\sqrt5\\quad\\textbf{(PROVED)}.$$
Conv-B (unit-normal) fold with $\\ell_j=\\log(A_j/|m_j|)=(\\log2,-\\frac12\\log2,0)$ — log 3 never appears —
equals $\\sqrt5\\,X$ with E5e's full 20-constant closed form $X$ **exactly** (all 20 coefficients), and matches the
independent 337-digit value $P^{(5)}_{arc}=-852.2552563361562115\\dots$ to 2.96e-123 (NUMERICAL confirmation).

## 4. Target 3 — the regular hexagon ($r=1$): PROVED

Paper's hexagon: $K_p=\\{\\sum_j|\\langle n_j,x\\rangle|^p\\le1\\}$ ($w_j=1$, edge $1/\\sqrt3$);
edge data $A_0=A_1=B=\\frac12,\\kappa=0,L_2=\\frac13,q_0=-\\frac16,V_{02}=\\frac13$, logs $\\log(A_j/w_j)=-\\log2$.
At $r=1$ all non-local terms vanish ($J\\sim x^c$); the local assembly (E5_polys × exact JLOC table) gives
$P^{(5)}(1)=6E_{arc}$, and **exactly**
$$\\frac{128}{\\sqrt3}P^{(5)}(1)=30720\\log2-4920\\pi^2-61440\\log^22+6240\\pi^2\\log2+30720\\log^32+448\\pi^4+5760\\zeta(3)\\log2-1920\\pi^2\\log^22-5120\\log^42+11340\\zeta(5)+5580\\zeta(2)\\zeta(3)-123\\pi^4\\log2-1440\\zeta(3)\\log^22+160\\pi^2\\log^32+256\\log^52,$$
with $[\\zeta(5)]P^{(5)}(1)=\\frac{2835\\sqrt3}{32}$ and $[\\zeta(2)\\zeta(3)]=\\frac{1395\\sqrt3}{32}$ exactly
(and no $\\log 3$ anywhere). Numeric: matches D17's independent 100-digit value to 4.98e-98.

## 5. Evidence table (all numbers printed by the listed code)

| Check | Object | Result | Residual | Digits | Code |
|---|---|---|---|---|---|
| relation rank vs DG | codims w1..5 | 1,2,3,5,8 | exact | exact | verify_rul1 / rul1_step02 |
| $\\zeta(\\bar2,1)$ | control | $\\zeta(3)/8$ | cert exact; 1.3e-121 | 120 | rul1_step03, verify |
| $\\zeta(\\bar3,\\bar1)$ | control | $-2Li_4+\\pi^4/180+\\pi^2\\log^22/12-\\log^42/12$ | cert exact; 1.3e-121 | 120 | rul1_step03, verify |
| five $J_{a,k}$ | Target 1 | D12 closed forms | certs exact; ≤ 9.7e-120 | 120 | rul1_step04, verify |
| $N_4,N_5$ closed forms | Target 1 | §2 formulas | exact vs D12 formula; 5.6e-132, 1.7e-131 | 130 | rul1_step04b, verify |
| 13 HPL atoms @$x=-1$ | Target 2 | §3 table | certs exact; ≤ 1.8e-120 | 120 | rul1_step05, verify |
| fold $\\zeta(5),\\zeta(2)\\zeta(3)$ | Target 2 | $-\\frac{20187}{2000}\\sqrt5$, $-\\frac{858}{125}\\sqrt5$ | **exact** | exact | rul1_step05b, verify |
| far-kink parts | Target 2 | $-\\frac{711}{400}\\sqrt5$, $-\\frac{1239}{250}\\sqrt5$ | **exact** | exact | rul1_step05b, verify |
| conv-B fold vs E5e $X$ | Target 2 | all 20 constants | exact; 2.96e-123 vs 337-digit data | 120+ | rul1_step05c, verify |
| regular hexagon | Target 3 | paper's 15-term form | exact; 4.98e-98 vs D17 100-digit | 98 | rul1_step06, verify |

## 6. Files

- `rul1.py` (engine: words/stuples, shuffle/stuffle, DS+DUAL relations) · `rul1la.py` (exact Gauss + history)
- `rul1_named.py` (named basis word-combos) · `rul1_reduce.py` (reductions + certificates) · `rul1num.py` (E1/E2/E3 evaluators)
- `rul1_step01.py` sanity · `rul1_step02.py` relation build + RREF + named inverses (`rul1_data.pkl`)
- `rul1_step03.py` controls · `rul1_step04.py` five J's (`rul1_certs_J.pkl`) · `rul1_step04b.py` $N_4,N_5$
- `rul1_step05.py` 13 atoms (`rul1_certs_atoms.pkl`) · `rul1_step05b.py` polar fold (ζ-coords)
- `rul1_step05c.py` conv-B fold vs E5e · `rul1_step05d.py` numeric convention check
- `rul1_step06.py` regular hexagon · **`verify_rul1.py` (41/41 PASS, 235 s)**
- `progress.log`

## 7. Label summary

- **PROVED:** the whole reduction system (DS+DUAL relations, exact ranks); the five J-values; $N_4,N_5$;
  the 13 atom values; the folded $\\zeta(5),\\zeta(2)\\zeta(3)$ coefficients $-\\frac{20187}{2000}\\sqrt5,-\\frac{858}{125}\\sqrt5$
  and far-kink parts $-\\frac{711}{400}\\sqrt5,-\\frac{1239}{250}\\sqrt5$; the conv-B arc value (= $\\sqrt5\\,X_{E5e}$, all 20 constants);
  the regular hexagon's 15-term form and its $\\frac{2835\\sqrt3}{32},\\frac{1395\\sqrt3}{32}$ coefficients.
  (In each case the proof = an exact rational identity modulo proved DS/DUAL relation instances, with a stored certificate.)
- **NUMERICAL:** every digit-agreement column (E1/E2/quadrature comparisons, E5e 337-digit data, D17 100-digit data);
  the linear independence of the basis constants (used only for "the" coordinates, certified by Deligne–Goncharov + exact dets).
'''
pathlib.Path('result.md').write_text(md)
print("result.md written", len(md.splitlines()), "lines")

import subprocess, time
log_line = time.strftime('%Y-%m-%d %H:%M:%S') + " RUL1 COMPLETE: verify_rul1.py 41/41 PASS (235s); Targets 1,2,3 PROVED"
with open('progress.log','a') as f:
    f.write(log_line + '\n')
print(subprocess.run(['bash','-c','ls -la; echo ===; tail -5 progress.log'], capture_output=True, text=True).stdout)
