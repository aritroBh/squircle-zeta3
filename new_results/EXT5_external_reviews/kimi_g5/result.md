# Referee report on G5 (Theorems 12 and 13, decay of $M_k$)

**Verdict: both theorems are correctly proved; the machine steps are reproducible and their outputs match the shipped `.out` files; the new main-term argument (Region I crude bound + Region II explicit saddle point) is sound.** I found three cosmetic numerical misstatements in the text and one inaccurate entry in the (non-load-bearing) checks table; none affects any proved bound. Details below.

---

## 1. What was checked

Files read from `<project>/new_results/EXT5_external_reviews/g5_pack/` (read-only): `main_tex_statements_l1255-1300.tex`, `sec_g5_decay.tex`, `sec_g5_varpi1.tex`, `macros_g5.tex`, `REPORT_G5.md`, `g5_consts.py/.out`, `g5_varpi1_algebra.py/.out`, `g5_checks.py/.out`, `wrap_g5.pdf`. All work was written only to my working directory (`steps/`).

**Statement match.** Theorem 12 (sec_g5_decay.tex): $|M_k+\frac{\sqrt\pi}{N}e^{-\theta}\cos\theta|\le 1.21\cdot10^9 N^{-3/2}e^{-\theta}$, $N=k+1\ge100$ (the user's shorthand "$|M_k-\mathrm{main}|$" = $M_k$ minus this cosine main term); Theorem 13 (sec_g5_varpi1.tex): $|M_k+\frac{\sqrt\pi}{N}e^{-\theta}\operatorname{Re}[e^{i\theta}(1+\varpi_1/\sqrt N)]|\le 1{,}000{,}668\,N^{-2}e^{-\theta}$, same range. The pack explicitly does *not* reprove the paper's stronger inputs ($254.23$; $1030\,N^{-1/2}e^{-\sqrt{2\pi N}}$) — it replaces them by its own Prop. main ($71{,}000$) and Prop. rest(b) ($120{,}000$). The closing arithmetic in the brief — $10{,}600{,}000\cdot114.096+7101.49<1.21\cdot10^9$ and $71{,}000+120{,}000=191{,}000\le1{,}000{,}668$ — is exactly what the scripts certify, and it is correct.

## 2. Machine steps (re-run, bit-compared)

| step | script | result |
|---|---|---|
| step_02_rerun_consts.py | g5_consts.py (Arb, 200 bits) | reproduces `g5_consts.out` exactly (modulo the `elapsed` line); all 17 internal `PROVED (Arb)` assertions pass |
| step_03_rerun_algebra.py | g5_varpi1_algebra.py (SymPy) | reproduces `g5_varpi1_algebra.out` exactly; exact identities asserted ($\nu_0$, $\eta_1$, $\varpi_1\sqrt{2/\pi}=\eta_1+(\nu_0-\tfrac18)\tfrac{1+i}{\sqrt2\pi}$, $\Lambda_N=-i\sqrt\pi/N$) |
| step_04_rerun_checks.py | g5_checks.py (mpmath, ~211 s) | reproduces `g5_checks.out` exactly (modulo timing) |

## 3. Independent recomputation (my own code, not the authors' scripts)

**step_05_verify_constants.py** (mpmath, 80 digits): every constant in `g5_consts.out` recomputed from the tex formulas — all agree to the printed enclosure digits, e.g. $B_1=70.09858533511162375$, $C_\star=10{,}578{,}907.823529397$, $(3/(e\delta))^3=114.096375458588$, $c_{\rm rest}=114{,}536.103195$, $S_1(100)=5.989464055$, $c_1=5.08862391583$, $E_{II}(10^4)=7912.190712134046$, Thm-12 total $1{,}209{,}428{,}681.354186$ (margin $571{,}318.6$ below $1.21\cdot10^9$). One nuance: the text's $c_I=70{,}217.69$ uses the relaxed $c_1\le5.1$; with the certified $c_1=5.0886$ it would be $70{,}103.93$ — both $<71{,}000$. Fine.

**step_06a_identities.py** (70 digits): exact identities verified numerically ($\nu_0,\eta_1,\varpi_1$-identity, $\Lambda_N$ at $N=10^2,10^4,10^7$; $2K_{3/2}(\theta)=\sqrt{2\pi/\theta}(1+1/\theta)e^{-\theta}$); grid checks — $\min|1+E^2|=0.7178\ge1/2$ on $\mathcal L$; $|v(u)-1+u/2|\le\tfrac32u^2$; Lemma H constant $\tfrac{45}8c^{-7/2}$ (worst grid ratio $0.032$); $(2m)!\ge24\cdot30^{m-2}$, $\ge720\cdot56^{m-3}$; $q(\sigma)$ bounds with the sharp constant $896/895$; $c_w=0.0591\ge0.0590486$, $c_{w'}=0.00197\ge0.0019664$. All OK.

**step_06b_integral_bounds.py** (quadrature): Lemma $F(j,q;\lambda)$ (worst ratio $0.074$), Lemma $W_k$ (worst ratio $0.315$), $\mu$-lemma (a) $0\le\operatorname{Re}R\le\gamma_0ye^{-3\sigma}$ and (c) tangent slopes — all OK.

**step_07a/07b/07c** (independent quadrature of $J$, $E_{II}$ and rest-piece monotonicity, the $g'$ slip): see findings below. Confirmed: $Ne^\theta|J|\le1.9387\ll5.0886$ for $N=100..10^4$; $N^2e^\theta|J-T|\le1.78\ll7912.2$; $E_{II}(N)$ nonincreasing on a log grid $10^4..10^{16}$ (e.g. $E_{II}=5687.9, 3423.0, 3288.1$ at $N=2\cdot10^4,10^5,10^{10}$); all six rest pieces $p_0..p_5$ nonincreasing for $N\in[100,10^7]$; $|g(N)-g(k)|\le1.9\,k^{-3/2}e^{-\sqrt{\pi k}}$ holds with worst ratio $0.663$ ($k=99..499$).

## 4. Hand audit of the proofs

I checked every lemma of both sections line by line (kernels lemma incl. $|1+E^2|\ge1/2$ via $|1+E^2|^2=(e^{-y}-\sin y)^2+\cos^2y$; logit form; kernel identity by contour indentation incl. the $V_J$ bound $1/(4J)$; $M_k$ representation with the $c_j$ majorant $m!(4\sqrt2/3)^m(\log2)^{j-m}$; rotation to $\mathcal L$; $W_k$ bound by per-term rotation; $F(j,q;\lambda)$ by AM–GM; rest assembly; Region I via $2K_{3/2}(\theta)=\sqrt{2\pi/\theta}(1+1/\theta)e^{-\theta}$ with $S_1$ monotonicity; $\mu$ convexity and secant slopes; $h(s)$ Watson bound with $\tfrac{45}8$; Steps 1–8 of Region II: factorisation $f=e^{-z}Ae^{-z(\cosh\sigma-1)}$, strip bounds ($E_1=0.07255<0.08$), Cauchy estimates ($a_4=256m$ correct via Taylor's integral form), window constants, Gaussian model, tails, exact identification of $\varpi_1$ and $\Lambda_N$, assembly and monotonicity). **No gap found.** Theorem 12: $C_\star(3/(e\delta))^3+7100+\sqrt\pi|\varpi_1|=1{,}209{,}428{,}681.35<1.21\cdot10^9$ ✓ (sup of $u^3e^{-\delta u}$ at $u=3/\delta=13.18>10$ ✓; exponent inequalities for part (a) need $N\ge18$ ✓). Theorem 13: $71{,}000+120{,}000=191{,}000\le1{,}000{,}668$ ✓.

## 5. Findings (all minor; none affects any proved statement)

1. **Cosmetic, proof of (Mk) in Thm 12:** the displayed "$|g'(t)|\le(t^{-2}+\sqrt\pi t^{-3/2})e^{-\sqrt{\pi t}}\le(0.1+\sqrt\pi)t^{-3/2}e^{-\sqrt{\pi t}}$ for $t\ge99$" is false for $t\in[99,100)$ ($t^{-1/2}\le0.1$ needs $t\ge100$; at $t=99$: coefficient $1/\sqrt{99}+\sqrt\pi=1.87296>0.1+\sqrt\pi=1.87245$, excess $5.1\cdot10^{-7}$, confirmed numerically). The final bound $|g(N)-g(k)|\le1.9\,k^{-3/2}e^{-\sqrt{\pi k}}$ that is actually used is correct (worst ratio $0.663$).
2. **Cosmetic, Region II Step 8 monotonicity:** the text says "$N^2\rho^{-3/2}/s_+\propto N^{3/4}e^{-\kappa\theta}$ and $N^2\rho^{-1}/(s_--1)\propto Ne^{-\kappa\theta}$"; the true powers are $N^{5/4}$ and $N^{3/2}$. With the correct powers the sign conditions are $\sqrt N>44.9$ and $\sqrt N>53.9$ (vs. the text's $26.9, 35.9$) — still $<100$, and $E_{II}$ is confirmed nonincreasing numerically up to $10^{16}$. Conclusion unaffected.
3. **Cosmetic, Step 4:** "$0\le q(\sigma)-\sigma^4/24\le\frac{1.0011}{720}\sigma^6$" — the sharp constant is $(897/896)/720=1.001116/720$, so the printed $1.0011$ is off by $1.5\cdot10^{-6}$ (relative). The constant actually used downstream, $c_{w'}=0.00197$, covers the true $\sqrt2(896/895)/720=0.0019664$. No effect.
4. **Numerical-evidence only:** `g5_checks.py`'s table entry $N^2e^\theta|J-T|=1.3095738$ at $N=10^4$ is quadrature-inaccurate; two independent finer quadratures of mine (60–80 digits, 400 and 1600 subintervals across the saddle window) agree with each other to 9 digits and give $1.5775346$ (at $N=3000$: authors $1.5957761$, mine $1.5931209$). All $N\le1000$ entries and all maxima over $N\in[100,420]$ reproduce correctly. The proofs use only Arb-certified constants plus hand-proved monotonicity, so this has no logical bearing.

## 6. Step files (in my working directory `steps/`)

`step_02_rerun_consts.py` (Arb rerun + diff), `step_03_rerun_algebra.py` (SymPy rerun + diff), `step_04_rerun_checks.py` (mpmath rerun + diff, 211 s), `step_05_verify_constants.py` (independent 80-digit recomputation of all constants), `step_06a_identities.py`, `step_06b_integral_bounds.py` (independent verification of all analytic inequalities), `step_07a_regions.py` (region bounds, $E_{II}$ and rest-piece monotonicity, $g'$ slip), `step_07b_Jfine.py`, `step_07c_Jcross.py` (fine quadrature establishing the true $N=3000,10^4$ check values; the 07d attempt at 60 dps with 3200 subintervals hit an mpmath quadrature instability and is superseded by 07b/07c, which agree with each other).